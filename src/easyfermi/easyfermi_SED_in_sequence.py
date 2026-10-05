"""
easyfermi_SED_in_sequence.py
============================

Helper module for easyfermi: computes the gamma-ray SED of the target for every time bin of a
pre-computed light curve (fermipy/easyfermi ``*_lightcurve.fits``), and gives access to the data of
all the ``*sed.fits`` files created (one per time bin).

No spectral fit is done here: the module only runs gta.sed() in each bin. The SED data are read back
from the fits files, so they can be fitted later (also in a different session, without recomputing).


Public functions
----------------
SED_in_sequence(...) -> runs gta.sed() in every time bin (optionally in parallel) and returns read_SEDs(...)
read_SEDs(...)       -> reads the *sed.fits files already existing in the time-bin directories
"""

import glob
import gc
import os
import platform
import re
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import yaml
from astropy.io import fits


# ----------------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------------
def _find_bin_directories(LC_directory):
    """Returns the time-bin directories (lightcurve_XXX) of a light curve, sorted in time."""
    directories = [d for d in glob.glob(os.path.join(LC_directory, "lightcurve*"))
                   if os.path.isdir(d) and os.path.isfile(os.path.join(d, "config.yaml"))]  # *.log files are skipped
    return sorted(directories, key=_bin_sort_key)


def _bin_name_numbers(bin_directory):
    """Integers at the end of the directory name: [tmin, tmax] (MET) for fermipy bins "lightcurve_<tmin>_<tmax>",
    or [n] for bins named "lightcurve_<n>"."""
    match = re.search(r"(\d+)(?:_(\d+))?$", os.path.basename(os.path.normpath(bin_directory)))
    if match is None:
        return []
    return [int(g) for g in match.groups() if g is not None]


def _bin_sort_key(bin_directory):
    numbers = _bin_name_numbers(bin_directory)
    return numbers[0] if numbers else -1


def _bin_number(bin_directory, tmin_met_LC=None):
    """
    Row of the light-curve table that corresponds to a time-bin directory.
    For fermipy bins "lightcurve_<tmin>_<tmax>", the row is found by matching tmin (MET) with the table;
    for bins "lightcurve_<n>", n is the row.
    """
    numbers = _bin_name_numbers(bin_directory)
    if len(numbers) == 2:
        if tmin_met_LC is None or len(tmin_met_LC) == 0:
            return -1
        row = int(np.argmin(np.abs(tmin_met_LC - numbers[0])))
        return row if abs(tmin_met_LC[row] - numbers[0]) < 2.0 else -1  # the directory names use integer seconds
    return numbers[0] if numbers else -1


def _sed_files(bin_directory):
    return sorted(glob.glob(os.path.join(bin_directory, "*sed.fits")))


def _as_float_array(column):
    return np.array(np.ma.filled(np.ma.asarray(column, dtype=float), np.nan), dtype=float)


def _table_to_dict(data):
    """Converts a FITS table into a dict {column_name: float array} (non-numeric columns are skipped)."""
    out = {}
    for name in data.names:
        try:
            out[name.lower()] = _as_float_array(data[name])
        except (TypeError, ValueError):
            pass
    return out


def _read_LC_times(LC_file):
    """Returns (tmin_mjd, tmax_mjd, tmin_met) of the light-curve bins (None where they cannot be read)."""
    try:
        with fits.open(LC_file) as hdul:
            names = [h.name.upper() for h in hdul]
            data = hdul["LIGHTCURVE"].data if "LIGHTCURVE" in names else hdul[1].data
            tmin_met = _as_float_array(data["tmin"]) if "tmin" in [n.lower() for n in data.names] else None
            return _as_float_array(data["tmin_mjd"]), _as_float_array(data["tmax_mjd"]), tmin_met
    except Exception:
        return None, None, None


# ----------------------------------------------------------------------------
# Reading the SEDs
# ----------------------------------------------------------------------------
def _read_one_SED(sed_file, bin_directory, tmin_LC, tmax_LC, tmin_met_LC=None):
    with fits.open(sed_file) as hdul:
        names = [h.name.upper() for h in hdul]
        sed = _table_to_dict(hdul["SED"].data if "SED" in names else hdul[1].data)
        model_flux = _table_to_dict(hdul["MODEL_FLUX"].data) if "MODEL_FLUX" in names else None

    if "e_ctr" not in sed and "e_min" in sed and "e_max" in sed:
        sed["e_ctr"] = np.sqrt(sed["e_min"] * sed["e_max"])

    index = _bin_number(bin_directory, tmin_met_LC)
    tmin = tmax = None
    if tmin_LC is not None and 0 <= index < len(tmin_LC):
        tmin, tmax = float(tmin_LC[index]), float(tmax_LC[index])

    return {
        "bin_index": index,
        "bin_directory": os.path.abspath(bin_directory),
        "sed_file": os.path.abspath(sed_file),
        "tmin_mjd": tmin,
        "tmax_mjd": tmax,
        "sed": sed,                  # e_min, e_max, e_ctr, e2dnde, e2dnde_err, e2dnde_ul95, ts, ... (MeV, MeV cm-2 s-1)
        "model_flux": model_flux,    # energies, dnde, dnde_lo, dnde_hi of the best-fit model (or None)
    }


def read_SEDs(which_LC):
    """
    Reads the *sed.fits files of all the time bins of a light curve.

    Parameters
    ----------
    which_LC : str
        Path to the light-curve fits file (e.g. .../Bayesian_blocks_light_curve_p0=0.05/bayesian_blocks_lightcurve.fits).
        The time-bin directories (lightcurve_XXX) must be in the same directory as this file.

    Returns
    -------
    SED_data : list of dict, sorted by time bin. Time bins without a *sed.fits file are not included.
        Each element has: bin_index, bin_directory, sed_file, tmin_mjd, tmax_mjd,
        sed (dict of arrays, one entry per column of the SED table) and model_flux (dict or None).
    """
    LC_file = os.path.abspath(which_LC)
    if not os.path.isfile(LC_file):
        raise FileNotFoundError(f"Light-curve file not found: {LC_file}")
    LC_directory = os.path.dirname(LC_file)
    tmin_LC, tmax_LC, tmin_met_LC = _read_LC_times(LC_file)

    SED_data = []
    for bin_directory in _find_bin_directories(LC_directory):
        files = _sed_files(bin_directory)
        if len(files) == 0:
            continue
        SED_data.append(_read_one_SED(files[0], bin_directory, tmin_LC, tmax_LC, tmin_met_LC))
    return SED_data


# ----------------------------------------------------------------------------
# SED of one time bin (this function runs in the worker processes)
# ----------------------------------------------------------------------------
def _SED_one_bin(bin_directory, target, emin, emax, number_of_bins, use_local_index, free_radius, overwrite):

    try:
        if not overwrite and len(_sed_files(bin_directory)) > 0:
            return bin_directory, "skipped", "SED already exists"

        from fermipy.gtanalysis import GTAnalysis  # imported here to keep the module light and fork-safe

        config_file = os.path.join(bin_directory, "config.yaml")
        with open(config_file, "r") as file:
            config = yaml.safe_load(file)
        config["fileio"]["outdir"] = str(bin_directory)
        config["fileio"]["logfile"] = str(bin_directory)  # fermipy adds the .log extension
        with open(config_file, "w") as file:
            yaml.dump(config, file)

        gta = GTAnalysis(config_file, logging={"verbosity": 3})
        gta.setup()

        if target is None:
            target = gta.config["selection"].get("target") or gta.roi.sources[0].name
        if emin is None:
            emin = gta.config["selection"]["emin"]
        if emax is None:
            emax = gta.config["selection"]["emax"]

        # Same free parameters used in the SED_in_Sequence.py script:
        gta.free_sources(distance=free_radius, pars="norm")
        for name in ("galdiff", "isodiff"):
            try:
                gta.free_source(name)
            except Exception:
                pass
        gta.free_source(target)

        fit_results = gta.fit()
        print(f"{os.path.basename(bin_directory)} - fit quality: {fit_results['fit_quality']}")

        loge_bins = np.linspace(np.log10(emin), np.log10(emax), num=number_of_bins)
        gta.sed(target, loge_bins=loge_bins, make_plots=False, use_local_index=use_local_index,
                write_fits=True, write_npy=False)

        del gta
        gc.collect()
        return bin_directory, "done", ""

    except Exception as error:
        return bin_directory, "failed", f"{type(error).__name__}: {error}"


# ----------------------------------------------------------------------------
# Main function
# ----------------------------------------------------------------------------
def SED_in_sequence(target, which_LC, number_of_bins=10, use_local_index=False, use_multiprocessing=True,
                    n_cores=None, emin=None, emax=None, free_radius=3.0,
                    overwrite=False, mp_start_method="fork"):
    """
    Computes the SED of the target in every time bin of a light curve.

    Parameters
    ----------
    target : str
        Name of the target source.
    which_LC : str
        Path to the light-curve fits file. The time-bin directories (lightcurve_XXX) must be in the
        same directory as this file.
    number_of_bins : int
        Number of energy bins of each SED (logarithmically spaced between emin and emax).
    use_local_index : bool
        If True, the spectral index used to compute the flux in each energy bin is the local index
        of the bin (gta.sed use_local_index). If False, the index of the global model is used.
    use_multiprocessing : bool
        If True, the time bins are processed in parallel (one process per time bin). Ignored (serial
        calculation) on macOS and Windows, where the "fork" start method is not safely available.
    n_cores : int or None
        Number of processes. If None, (number of CPUs - 1) is used.
    emin, emax : float or None
        Energy range of the SEDs in MeV. If None, the range of the analysis (config.yaml) is used.
    free_radius : float
        Sources within this radius (deg) from the ROI center have their normalization freed in the fit.
    overwrite : bool
        If False, time bins that already have a *sed.fits file are skipped.
    mp_start_method : str
        Multiprocessing start method ("fork" on Linux).

    Returns
    -------
    SED_data : list of dict
        Output of read_SEDs(which_LC).
    """
    LC_file = os.path.abspath(which_LC)
    if not os.path.isfile(LC_file):
        raise FileNotFoundError(f"Light-curve file not found: {LC_file}")

    bin_directories = _find_bin_directories(os.path.dirname(LC_file))
    if len(bin_directories) == 0:
        raise FileNotFoundError(f"No lightcurve_XXX directory (with config.yaml) found in {os.path.dirname(LC_file)}.")

    kwargs = dict(target=target, emin=emin, emax=emax, number_of_bins=number_of_bins,
                  use_local_index=use_local_index, free_radius=free_radius, overwrite=overwrite)

    if use_multiprocessing and platform.system() != "Linux":
        print(f"Multiprocessing is only used on Linux (this system: {platform.system()}). Running without it.")
        use_multiprocessing = False

    print(f"Computing the SEDs of {len(bin_directories)} time bins...")
    results = []
    if use_multiprocessing and len(bin_directories) > 1:
        if n_cores is None:
            n_cores = max(1, (os.cpu_count() or 2) - 1)
        n_cores = int(max(1, min(n_cores, len(bin_directories))))
        with ProcessPoolExecutor(max_workers=n_cores, mp_context=mp.get_context(mp_start_method)) as executor:
            futures = [executor.submit(_SED_one_bin, d, **kwargs) for d in bin_directories]
            for future in as_completed(futures):
                results.append(future.result())
                print(f"  {os.path.basename(results[-1][0])}: {results[-1][1]} {results[-1][2]}")
    else:
        for d in bin_directories:
            results.append(_SED_one_bin(d, **kwargs))
            print(f"  {os.path.basename(results[-1][0])}: {results[-1][1]} {results[-1][2]}")

    n_done = sum(r[1] == "done" for r in results)
    n_skipped = sum(r[1] == "skipped" for r in results)
    failed = [os.path.basename(r[0]) for r in results if r[1] == "failed"]
    print(f"SEDs: {n_done} computed, {n_skipped} skipped, {len(failed)} failed" + (f" ({', '.join(failed)})." if failed else "."))

    return read_SEDs(LC_file)


# ----------------------------------------------------------------------------
# Quick command-line test:  python easyfermi_SED_in_sequence.py path/to/lightcurve.fits [number_of_bins]
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        sys.exit("Usage: python easyfermi_SED_in_sequence.py <lightcurve.fits> [number_of_bins]")
    _N = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    _data = SED_in_sequence(sys.argv[1], number_of_bins=_N)
    for _d in _data:
        print(_d["bin_index"], _d["tmin_mjd"], _d["tmax_mjd"], _d["sed_file"])
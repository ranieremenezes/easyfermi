"""
easyfermi_bayesian_blocks.py
============================

Helper module for easyfermi: reads a gamma-ray light curve (.csv or .fits) and
computes its Bayesian-blocks representation (Scargle et al. 2013, ApJ 764, 167),
using the "point measurements" fitness function (flux +- error in each time bin).

Supported input formats
-----------------------
* FITS: the light-curve file written by fermipy / easyfermi
  (``TARGET_NAME_lightcurve.fits``, table extension "LIGHTCURVE").
  Columns used: tmin_mjd, tmax_mjd, ts, flux/eflux, flux_err/eflux_err,
  flux_ul95/eflux_ul95.

* CSV: the file downloaded from the Fermi-LAT Light Curve Repository (LCR).
  Columns used: MET (or Julian Date), TS, Photon Flux, Photon Flux Error.
  In this format, non-detections come as "< value" (upper limit) with no
  flux/error, and the bin width is not given, so it is inferred from the cadence
  (or set by the user through ``bin_width_days``).

Public functions
----------------
read_light_curve(...)            -> dict with the light-curve data
restrict_light_curve(...)        -> keeps only the bins inside the time range of the data
bayesian_blocks_light_curve(...) -> dict with the Bayesian blocks
MET_to_MJD(...), MJD_to_MET(...) -> time conversions (fermipy convention)
"""

import os

import numpy as np
from astropy.io import fits
from astropy.stats import bayesian_blocks
from astropy.table import Table

# MJD (TT) corresponding to Fermi MET = 0 (same convention used by fermipy)
MJDREF_MET = 51910.0007428703703703703
SECONDS_PER_DAY = 86400.0

_FLUX_UNITS = {
    "flux": "ph cm-2 s-1",
    "eflux": "MeV cm-2 s-1",
}


# ----------------------------------------------------------------------------
# Reading
# ----------------------------------------------------------------------------
def _guess_format(input_file):
    ext = os.path.splitext(str(input_file))[1].lower()
    if ext == ".csv":
        return "csv"
    if ext in (".fits", ".fit", ".fts"):
        return "fits"
    raise ValueError(
        f"Could not infer the file format from the extension '{ext}'. "
        "Use file_format='csv' or file_format='fits'."
    )


def _as_float_array(column):
    """Converts a (possibly big-endian / masked) FITS or astropy column to a native float array."""
    return np.array(np.ma.filled(np.ma.asarray(column, dtype=float), np.nan), dtype=float)


def _find_column(names, *prefixes):
    """Returns the first column whose (lower-case, stripped) name starts with one of the prefixes."""
    lowered = {name.strip().lower(): name for name in names}
    for prefix in prefixes:
        for low, original in lowered.items():
            if low.startswith(prefix.lower()):
                return original
    raise KeyError(f"None of the columns {prefixes} were found. Available columns: {list(names)}")


def _parse_lcr_column(values):
    """
    Parses a Light Curve Repository flux column.

    Entries look like '5.85e-8' (measurement), '< 5.83e-8' (upper limit) or '-' (no value).

    Returns
    -------
    measurement, upper_limit : arrays of floats (NaN where not applicable)
    """
    values = np.asarray(values).astype(str)
    measurement = np.full(len(values), np.nan)
    upper_limit = np.full(len(values), np.nan)
    for i, raw in enumerate(values):
        text = raw.strip().strip('"')
        try:
            if text.startswith("<"):
                upper_limit[i] = float(text[1:])
            else:
                measurement[i] = float(text)
        except ValueError:  # '-', '', 'nan', ...
            pass
    return measurement, upper_limit


def _read_fits(input_file, flux_type):
    with fits.open(input_file) as hdul:
        if "LIGHTCURVE" in [h.name.upper() for h in hdul]:
            data = hdul["LIGHTCURVE"].data
        else:
            data = hdul[1].data

        needed = ["tmin_mjd", "tmax_mjd", "ts", flux_type, flux_type + "_err", flux_type + "_ul95"]
        missing = [c for c in needed if c not in data.names]
        if missing:
            raise KeyError(f"Columns {missing} not found in {input_file}.")

        out = {
            "tmin_mjd": _as_float_array(data["tmin_mjd"]),
            "tmax_mjd": _as_float_array(data["tmax_mjd"]),
            "ts": _as_float_array(data["ts"]),
            "flux": _as_float_array(data[flux_type]),
            "flux_err": _as_float_array(data[flux_type + "_err"]),
            "flux_ul95": _as_float_array(data[flux_type + "_ul95"]),
        }
    return out


def _read_csv(input_file, flux_type, bin_width_days):
    if flux_type != "flux":
        raise ValueError(
            "Light Curve Repository .csv files only contain photon flux. "
            "Use flux_type='flux' for this file."
        )

    table = Table.read(input_file, format="ascii.csv")
    names = table.colnames

    col_ts = _find_column(names, "ts")
    col_flux = _find_column(names, "photon flux [", "photon flux")
    col_err = _find_column(names, "photon flux error")

    ts = _as_float_array(table[col_ts])
    flux, flux_ul = _parse_lcr_column(table[col_flux])
    flux_err, _ = _parse_lcr_column(table[col_err])

    # Start of each bin: MET is preferred (exact); Julian Date is only rounded to 1 day in the LCR files.
    try:
        col_met = _find_column(names, "met")
        tmin_mjd = MJDREF_MET + _as_float_array(table[col_met]) / SECONDS_PER_DAY
    except KeyError:
        col_jd = _find_column(names, "julian date")
        tmin_mjd = _as_float_array(table[col_jd]) - 2400000.5

    order = np.argsort(tmin_mjd)
    tmin_mjd, ts, flux, flux_err, flux_ul = (a[order] for a in (tmin_mjd, ts, flux, flux_err, flux_ul))

    # Bin width: the LCR file does not give it. Some bins can be missing (failed fits), so the
    # nominal cadence is taken as the median spacing between consecutive bins.
    if bin_width_days is None:
        if len(tmin_mjd) < 2:
            raise ValueError("Only one time bin in the file: please provide bin_width_days.")
        bin_width_days = float(np.median(np.diff(tmin_mjd)))

    return {
        "tmin_mjd": tmin_mjd,
        "tmax_mjd": tmin_mjd + bin_width_days,
        "ts": ts,
        "flux": flux,
        "flux_err": flux_err,
        "flux_ul95": flux_ul,
    }


def read_light_curve(input_file, file_format="auto", flux_type="flux", bin_width_days=None):
    """
    Reads a light curve from a .csv (Fermi-LAT Light Curve Repository) or .fits (fermipy/easyfermi) file.

    Parameters
    ----------
    input_file : str
        Path to the light-curve file.
    file_format : {"auto", "csv", "fits"}
        File format. "auto" uses the file extension.
    flux_type : {"flux", "eflux"}
        "flux": photon flux [ph cm-2 s-1]; "eflux": energy flux [MeV cm-2 s-1] (only for .fits).
    bin_width_days : float or None
        Only used for .csv. Time-bin width in days. If None, it is inferred from the cadence.

    Returns
    -------
    lc : dict
        tmin_mjd, tmax_mjd, tmid_mjd, ts, flux, flux_err, flux_ul95 : arrays sorted in time.
        Values that do not exist are NaN (e.g. in .csv files, flux and flux_err are NaN for
        upper limits, and flux_ul95 is NaN for detections).
        Also: flux_type, flux_unit, file_format, input_file.
    """
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Light-curve file not found: {input_file}")

    file_format = str(file_format).lower()
    if file_format == "auto":
        file_format = _guess_format(input_file)
    if file_format not in ("csv", "fits"):
        raise ValueError("file_format must be 'auto', 'csv' or 'fits'.")

    flux_type = str(flux_type).lower()
    if flux_type not in _FLUX_UNITS:
        raise ValueError("flux_type must be 'flux' or 'eflux'.")

    if file_format == "fits":
        lc = _read_fits(input_file, flux_type)
        order = np.argsort(lc["tmin_mjd"])
        lc = {k: v[order] for k, v in lc.items()}
    else:
        lc = _read_csv(input_file, flux_type, bin_width_days)

    lc["tmid_mjd"] = 0.5 * (lc["tmin_mjd"] + lc["tmax_mjd"])
    lc["flux_type"] = flux_type
    lc["flux_unit"] = _FLUX_UNITS[flux_type]
    lc["file_format"] = file_format
    lc["input_file"] = os.path.abspath(input_file)
    return lc


# ----------------------------------------------------------------------------
# Time range of the data
# ----------------------------------------------------------------------------
def MET_to_MJD(met):
    """Fermi Mission Elapsed Time (s) -> MJD (same convention used by fermipy)."""
    return MJDREF_MET + np.asarray(met, dtype=float) / SECONDS_PER_DAY


def MJD_to_MET(mjd):
    """MJD -> Fermi Mission Elapsed Time (s) (same convention used by fermipy)."""
    return (np.asarray(mjd, dtype=float) - MJDREF_MET) * SECONDS_PER_DAY


_PER_BIN_KEYS = ("tmin_mjd", "tmax_mjd", "tmid_mjd", "ts", "flux", "flux_err", "flux_ul95")


def restrict_light_curve(lc, tmin_mjd, tmax_mjd):
    """
    Keeps only the time bins of a light curve that overlap the interval [tmin_mjd, tmax_mjd]
    (e.g. the time range of the Fermi-LAT data being analyzed).

    A bin that is only partially inside the interval is kept: its measurement still carries
    information about the variability inside the interval. The edges of the final Bayesian-blocks
    light curve must therefore be clipped to [tmin, tmax] afterwards.

    Parameters
    ----------
    lc : dict
        Output of read_light_curve().
    tmin_mjd, tmax_mjd : float
        Time interval (MJD) where data are available.

    Returns
    -------
    lc_restricted : dict
        Same structure as ``lc``, with only the overlapping bins. It also contains data_tmin_mjd,
        data_tmax_mjd, n_bins_original and n_bins_removed.
    inside : bool array (same length as the arrays of ``lc``)
        True for the bins kept.

    Raises
    ------
    ValueError
        If no bin of the light curve overlaps the interval.
    """
    tmin_mjd, tmax_mjd = float(tmin_mjd), float(tmax_mjd)
    if tmax_mjd <= tmin_mjd:
        raise ValueError(f"Invalid time interval: tmin_mjd = {tmin_mjd} >= tmax_mjd = {tmax_mjd}.")

    inside = (lc["tmax_mjd"] > tmin_mjd) & (lc["tmin_mjd"] < tmax_mjd)
    if not inside.any():
        raise ValueError(
            f"The light curve (MJD {np.min(lc['tmin_mjd']):.2f} - {np.max(lc['tmax_mjd']):.2f}) does not overlap "
            f"the time range of the data (MJD {tmin_mjd:.2f} - {tmax_mjd:.2f})."
        )

    lc_restricted = {key: (value[inside] if key in _PER_BIN_KEYS else value) for key, value in lc.items()}
    lc_restricted["data_tmin_mjd"] = tmin_mjd
    lc_restricted["data_tmax_mjd"] = tmax_mjd
    lc_restricted["n_bins_original"] = int(len(inside))
    lc_restricted["n_bins_removed"] = int((~inside).sum())
    return lc_restricted, inside


# ----------------------------------------------------------------------------
# Bayesian blocks
# ----------------------------------------------------------------------------
def bayesian_blocks_light_curve(lc, ts_min=4.0, p0=0.05, ncp_prior=None):
    """
    Computes the Bayesian blocks of a light curve read with read_light_curve().

    Only time bins with a valid measurement (finite flux and error > 0) and TS >= ts_min enter
    the calculation; the remaining bins (upper limits) are ignored, but they stay in ``lc``.
    A block therefore starts at the beginning of its first bin and ends at the end of its last bin,
    and it can span ignored bins located between them.

    Parameters
    ----------
    lc : dict
        Output of read_light_curve().
    ts_min : float or None
        Minimum TS for a bin to be used. If None, every bin with valid flux and error is used
        (in .csv files, bins with upper limits never have flux/error, so they are always skipped).
    p0 : float
        False-alarm probability of the Bayesian blocks (Scargle et al. 2013). Lower p0 -> fewer blocks.
    ncp_prior : float or None
        Prior on the number of change points. If given, it overrides p0.

    Returns
    -------
    bb : dict
        used                : bool array (same length as lc arrays), True for the bins used.
        block_id            : int array (bins used), block that each used bin belongs to.
        n_used, n_blocks    : number of bins used and number of blocks.
        block_tstart_mjd, block_tstop_mjd, block_tmid_mjd, block_duration_days : arrays (n_blocks).
        block_flux, block_flux_err : weighted mean and its error in each block (same unit as lc["flux"]).
        block_nbins         : number of bins used in each block.
        block_chi2_red      : reduced chi2 of each block against a constant (NaN for 1-bin blocks).
        flux_block_per_bin  : block flux assigned to each used bin.
        ts_min, p0, ncp_prior : parameters used.
    """
    flux = lc["flux"]
    err = lc["flux_err"]
    ts = lc["ts"]

    used = np.isfinite(flux) & np.isfinite(err) & (err > 0)
    if ts_min is not None:
        used &= np.isfinite(ts) & (ts >= ts_min)

    n_used = int(used.sum())
    if n_used == 0:
        raise ValueError("No valid time bin to compute the Bayesian blocks (check ts_min and the input file).")

    t = lc["tmid_mjd"][used]
    tmin = lc["tmin_mjd"][used]
    tmax = lc["tmax_mjd"][used]
    x = flux[used]
    sigma = err[used]

    if n_used == 1:
        block_id = np.zeros(1, dtype=int)
    else:
        # The fitness is scale invariant; rescaling only avoids huge 1/sigma^2 numbers (~1e16).
        scale = np.median(np.abs(x))
        scale = scale if scale > 0 else 1.0
        kwargs = {"ncp_prior": ncp_prior} if ncp_prior is not None else {"p0": p0}
        edges = bayesian_blocks(t, x / scale, sigma / scale, fitness="measures", **kwargs)
        # edges[1:-1] are the change points (mid-way between consecutive bins)
        block_id = np.searchsorted(edges[1:-1], t)

    n_blocks = int(block_id.max()) + 1

    b_tstart = np.empty(n_blocks)
    b_tstop = np.empty(n_blocks)
    b_flux = np.empty(n_blocks)
    b_err = np.empty(n_blocks)
    b_nbins = np.empty(n_blocks, dtype=int)
    b_chi2 = np.full(n_blocks, np.nan)
    flux_per_bin = np.empty(n_used)

    for k in range(n_blocks):
        m = block_id == k
        w = 1.0 / sigma[m] ** 2
        mean = np.sum(w * x[m]) / np.sum(w)
        b_tstart[k] = tmin[m].min()
        b_tstop[k] = tmax[m].max()
        b_flux[k] = mean
        b_err[k] = 1.0 / np.sqrt(np.sum(w))
        b_nbins[k] = int(m.sum())
        if b_nbins[k] > 1:
            b_chi2[k] = np.sum(w * (x[m] - mean) ** 2) / (b_nbins[k] - 1)
        flux_per_bin[m] = mean

    return {
        "used": used,
        "block_id": block_id,
        "n_used": n_used,
        "n_blocks": n_blocks,
        "block_tstart_mjd": b_tstart,
        "block_tstop_mjd": b_tstop,
        "block_tmid_mjd": 0.5 * (b_tstart + b_tstop),
        "block_duration_days": b_tstop - b_tstart,
        "block_flux": b_flux,
        "block_flux_err": b_err,
        "block_nbins": b_nbins,
        "block_chi2_red": b_chi2,
        "flux_block_per_bin": flux_per_bin,
        "ts_min": ts_min,
        "p0": p0,
        "ncp_prior": ncp_prior,
    }


def block_edges_mjd(lc, bb):
    """
    Edges (MJD, n_blocks + 1 values) of the time bins defined by the Bayesian blocks.

    - Internal edges: the edge between blocks k and k+1 is placed in the middle of the gap between the
      last bin used in block k and the first bin used in block k+1 (Scargle et al. 2013). The bins inside
      the gap (e.g. upper limits, not used in the calculation) are therefore split half to the block on the
      left and half to the block on the right. Without a gap, the edge is exactly the block boundary.
    - Outer edges: the limits of the light curve, so leading/trailing upper limits are absorbed by the
      first/last blocks.

    The edges are not clipped: clip them to the time range of the data if needed.
    """
    start, stop = bb["block_tstart_mjd"], bb["block_tstop_mjd"]
    return np.concatenate([[np.min(lc["tmin_mjd"])], 0.5 * (stop[:-1] + start[1:]), [np.max(lc["tmax_mjd"])]])


# ----------------------------------------------------------------------------
# Quick command-line test:  python easyfermi_bayesian_blocks.py file.csv|file.fits [ts_min] [p0]
# ----------------------------------------------------------------------------
if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        sys.exit("Usage: python easyfermi_bayesian_blocks.py <file.csv|file.fits> [ts_min] [p0]")
    _ts_min = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    _p0 = float(sys.argv[3]) if len(sys.argv) > 3 else 0.05
    _lc = read_light_curve(sys.argv[1])
    _bb = bayesian_blocks_light_curve(_lc, ts_min=_ts_min, p0=_p0)
    print(f"{len(_lc['ts'])} bins read, {_bb['n_used']} used, {_bb['n_blocks']} Bayesian blocks")
    for _k in range(_bb["n_blocks"]):
        print(
            f"  block {_k + 1:3d}: MJD {_bb['block_tstart_mjd'][_k]:.2f} - {_bb['block_tstop_mjd'][_k]:.2f}"
            f"  flux = {_bb['block_flux'][_k]:.3e} +- {_bb['block_flux_err'][_k]:.3e}"
            f"  ({_bb['block_nbins'][_k]} bins)"
        )
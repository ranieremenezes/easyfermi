"""
easyfermi_MCMC.py
=================

Helper module for easyfermi: MCMC (emcee) fit of gamma-ray SEDs, optional EBL correction (gammapy),
output of the results to the *sed.fits files, and quick-look plots.

It is used by:
    - Ui_mainWindow.EBL_and_MCMC()  -> MCMC fit of the SED of the full time interval (+ optional VHE data);
    - fit_SEDs_per_bin()            -> MCMC fit of the SED of every time bin of a light curve
                                       (data produced by Ui_mainWindow.SED_per_LC_bin()).

Conventions
-----------
- Energies in MeV, dN/dE in MeV-1 cm-2 s-1, E^2 dN/dE in MeV cm-2 s-1.
- Likelihood (SEDLikelihood): chi^2 in log10(dN/dE) for the data points (TS > TSmin) + the energy bins without
  detection (upper limits), as gaussian terms in linear flux built from their best-fit value and error (or as
  censored terms when only the upper limit is known). The fermipy likelihood scans (norm_scan/dloglike_scan) are
  not used. The old fit, with the upper limits ignored, is available with method="chi2".
- MCMC (run_MCMC): adaptive emcee run, vectorized over the walkers, that stops when the chain is longer than
  50 autocorrelation times; burn-in and thinning are set from the autocorrelation time.
- Priors: uniform within open intervals.

Public API
----------
MODELS, get_model(name)
SEDData                       -> container of the data points / upper limits (+ EBL-corrected versions)
EBL_absorption_model(...)     -> gammapy EBL absorption model
SEDLikelihood(...)            -> log-likelihood of a model given a SED (chi^2 + upper limits, or chi^2 only)
run_MCMC(likelihood) -> MCMCResult -> adaptive emcee run; best-fit parameters, posteriors and diagnostics
read_VHE_SED(...)             -> reads a VHE SED (gammapy format)
EBL_corrected_SED_columns(...), VHE_EBL_table(...), update_SED_fits(...)
corner_plot(...), plot_SED_MCMC(...)
fit_SEDs_per_bin(...)         -> MCMC of every SED of a light curve
"""

import os
import re
from dataclasses import dataclass
from typing import Callable, Optional, Sequence, Tuple

import numpy as np
import astropy.io.fits as pyfits
from astropy import units as u
import emcee
import corner
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator


# ----------------------------------------------------------------------------
# Default settings (same values used historically by easyfermi)
# ----------------------------------------------------------------------------
TS_MIN = 9.0            # SED bins with TS > TS_MIN are data points, the others are upper limits
QUANTILES = (0.16, 0.5, 0.84)

# Adaptive MCMC (see run_MCMC):
N_WALKERS = 50          # number of walkers (at least 2 x number of parameters + 2)
CHUNK_STEPS = 500       # the autocorrelation time is re-estimated after each chunk of steps
N_TAU_STOP = 50         # convergence: chain longer than N_TAU_STOP autocorrelation times...
TAU_RTOL = 0.01         # ... and autocorrelation time changing by less than 1% between chunks (emcee docs)
BURN_IN_TAU = 5         # burn-in = BURN_IN_TAU autocorrelation times
MAX_STEPS = 30000       # maximum number of steps per walker (the result is flagged if not converged)
INITIAL_SPREAD = 1e-3   # initial ball around the maximum-likelihood point, as a fraction of the prior width
N_PRIOR_SAMPLES = 4000  # random points of the prior used to start the maximum-likelihood search
MIN_N_TAU = 20          # chains shorter than MIN_N_TAU autocorrelation times are flagged as unreliable
UL_YERR_FRACTION = 0.3  # length of the upper-limit arrows, as a fraction of the UL

LOG10_E = np.log10(np.e)

# Labels of the GUI combo box -> name of the gammapy built-in EBL model
EBL_MODELS = {
    "Dominguez et al. (2011)": "dominguez",
    "Franceschini et al. (2008)": "franceschini",
    "Franceschini & Rodighiero (2017)": "franceschini17",
    "Saldana-Lopez et al. (2021)": "saldana-lopez21",
    "Finke et al. (2010)": "finke",
    "Finke et al. (2022) model A": "finke2022",
}


# ----------------------------------------------------------------------------
# Spectral models (all in log10 space; x = log10(E/MeV), log_e0 = log10(Emin/MeV) is the pivot energy)
# ----------------------------------------------------------------------------
def _power_law(theta, x, log_e0):
    N0, alpha = theta
    return N0 - alpha * (x - log_e0)


def _log_parabola(theta, x, log_e0):
    # dN/dE = N0 (E/E0)^(-alpha - beta ln(E/E0)), as in the Fermi-LAT LogParabola
    N0, alpha, beta = theta
    return N0 + (-alpha - beta * np.log(10.0) * (x - log_e0)) * (x - log_e0)


def _log_parabola_MTT(theta, x, log_e0):
    # log10(E^2 dN/dE) = log10(Sp) - alpha * log10(E/Ep)^2
    Splog, alpha, Ep = theta
    return Splog - alpha * (x - Ep) ** 2


def _exp_cutoff(x, Ec, b):
    # log10(exp(-(E/Ec)^b)), written without np.exp to avoid log10(0) = -inf at very high energies
    return -LOG10_E * (10.0 ** (x - Ec)) ** b


def _PLEC(theta, x, log_e0):
    N0, alpha, Ec, b = theta
    return N0 - alpha * (x - log_e0) + _exp_cutoff(x, Ec, b)


def _PLEC_bfix(theta, x, log_e0):
    N0, alpha, Ec = theta
    return N0 - alpha * (x - log_e0) + _exp_cutoff(x, Ec, 1.0)


def _PLEC_deMenezes(theta, x, log_e0):
    # log10(E^2 dN/dE), with Sp = peak flux, Ep = peak energy
    Sp, alpha, Ep, b = theta
    return Sp + (alpha - 2) * (Ep - x) + LOG10_E * ((2 - alpha) / b) * (1 - (10.0 ** (x - Ep)) ** b)


@dataclass(frozen=True)
class SpectralModel:
    """Definition of a spectral model fitted with the MCMC."""
    name: str
    function: Callable                           # f(theta, x, log_e0) -> log10(dN/dE) or log10(E^2 dN/dE)
    initial: Tuple[float, ...]                    # center of the initial gaussian ball of the walkers
    bounds: Tuple[Tuple[float, float], ...]       # uniform priors (open intervals)
    parameter_names: Tuple[str, ...]              # used in Target_results.txt and in the "MCMC Parameters" table
    posterior_names: Tuple[str, ...]              # columns of the "MCMC Posterior dist." table
    corner_labels: Tuple[str, ...]
    fixed_parameters: Tuple[Tuple[str, float], ...] = ()  # (name, value) of the fixed parameters
    uses_pivot: bool = True                       # if True, the row "Ep=Emin (log scale)" is saved
    fits_e2dnde: bool = False                     # True: function returns log10(E^2 dN/dE); False: log10(dN/dE)

    @property
    def ndim(self):
        return len(self.initial)

    def log_prior(self, theta):
        for value, (low, high) in zip(theta, self.bounds):
            if not low < value < high:  # also rejects NaN
                return -np.inf
        return 0.0

    def log_e2dnde(self, theta, x, log_e0):
        """log10(E^2 dN/dE) of the model, for any of the two parametrizations."""
        y = self.function(theta, x, log_e0)
        return y if self.fits_e2dnde else y + 2 * x


MODELS = {
    "PowerLaw": SpectralModel(
        "PowerLaw", _power_law, initial=(-13, 2.0), bounds=((-15, -7), (0.5, 5.0)),
        parameter_names=("N0 (log scale)", "Alpha"),
        posterior_names=("N0 distribution (log scale)", "Alpha distribution"),
        corner_labels=("N0", "alpha")),
    "LogPar": SpectralModel(
        "LogPar", _log_parabola, initial=(-13, 1.7, 0.2), bounds=((-15, -7), (1.0, 4.0), (-1, 1.0)),
        parameter_names=("N0 (log scale)", "Alpha", "Beta"),
        posterior_names=("N0 distribution (log scale)", "Alpha distribution", "Beta distribution"),
        corner_labels=("N0", "alpha", "beta")),
    "LogPar_MTT": SpectralModel(
        "LogPar_MTT", _log_parabola_MTT, initial=(-4.5, 0.2, 3.5), bounds=((-7, -1), (-1.0, 1.0), (2, 7)),
        parameter_names=("Sp (log scale)", "Alpha", "Ep (log scale)"),
        posterior_names=("Sp distribution (log scale)", "Alpha distribution", "Ep distribution (log scale)"),
        corner_labels=("Sp_log", "alpha", "Ep_log"), uses_pivot=False, fits_e2dnde=True),
    "PLEC": SpectralModel(
        "PLEC", _PLEC, initial=(-13, 1.7, 5, 1), bounds=((-15, -7), (1.0, 4.0), (3.0, 7.0), (0.2, 3.0)),
        parameter_names=("N0 (log scale)", "Alpha", "Ec", "b"),
        posterior_names=("N0 distribution (log scale)", "Alpha distribution", "Ec distribution", "b distribution"),
        corner_labels=("N0", "alpha", "Ec", "b")),
    "PLEC_bfix": SpectralModel(
        "PLEC_bfix", _PLEC_bfix, initial=(-13, 1.7, 5), bounds=((-15, -7), (1.0, 4.0), (3.0, 7.0)),
        parameter_names=("N0 (log scale)", "Alpha", "Ec"),
        posterior_names=("N0 distribution (log scale)", "Alpha distribution", "Ec distribution"),
        corner_labels=("N0", "alpha", "Ec"), fixed_parameters=(("b", 1.0),)),
    "PLEC_deMenezes": SpectralModel(
        "PLEC_deMenezes", _PLEC_deMenezes, initial=(-4, 1.7, 5, 1),
        bounds=((-8, -1), (0, 4.0), (2.0, 7.0), (0.01, 3.0)),
        parameter_names=("Sp (log scale)", "Alpha", "Ep (log scale)", "b"),
        posterior_names=("Sp distribution (log scale)", "Alpha distribution", "Ep distribution (log scale)", "b distribution"),
        corner_labels=("Sp_log", "alpha", "Ep", "b"), uses_pivot=False, fits_e2dnde=True),
}


def get_model(name):
    try:
        return MODELS[name]
    except KeyError:
        raise ValueError(f"Unknown MCMC model '{name}'. Available models: {', '.join(MODELS)}.") from None


# ----------------------------------------------------------------------------
# EBL
# ----------------------------------------------------------------------------
def EBL_absorption_model(EBL_label, redshift, EBL_path):
    """
    Returns the gammapy EBLAbsorptionNormSpectralModel for a given GUI label (see EBL_MODELS) and redshift.
    EBL_path is the directory containing the "ebl" folder with the EBL tables (easyfermi/resources).
    """
    from gammapy.modeling.models import EBLAbsorptionNormSpectralModel, EBL_DATA_BUILTIN  # heavy import

    try:
        EBL_model = EBL_MODELS[EBL_label]
    except KeyError:
        raise ValueError(f"Unknown EBL model '{EBL_label}'. Available models: {', '.join(EBL_MODELS)}.") from None
    if EBL_model == "finke2022":
        EBL_DATA_BUILTIN["finke2022"] = "$GAMMAPY_DATA/ebl/ebl_finke22.fits.gz"

    os.environ["GAMMAPY_DATA"] = str(EBL_path)  # gammapy will look for EBL models in this directory
    return EBLAbsorptionNormSpectralModel.read_builtin(EBL_model, redshift=redshift)


def EBL_attenuation(absorption, energy, redshift):
    """Attenuation factor exp(-tau) at the energies (MeV) given, as a float array."""
    energy = np.atleast_1d(np.asarray(energy, dtype=float))
    if energy.size == 0:
        return np.zeros(0)
    value = absorption.evaluate(energy * u.MeV, redshift, alpha_norm=1)
    return np.asarray(u.Quantity(value).to_value(u.dimensionless_unscaled), dtype=float)


# ----------------------------------------------------------------------------
# Data container
# ----------------------------------------------------------------------------
def _as_float(values):
    return np.asarray(values, dtype=float)


def _get_column(table, *names):
    """Returns the first column found (case insensitive) among `names` in a dict-like table, or None."""
    lower = {str(key).lower(): key for key in table.keys()}
    for name in names:
        if name.lower() in lower:
            return _as_float(table[lower[name.lower()]])
    return None


@dataclass
class SEDData:
    """
    Data points (TS > TSmin) and upper limits of a SED.
    The last n_VHE points (and n_VHE_UL upper limits) can be from a VHE instrument.
    """
    energy: np.ndarray                  # MeV
    energy_err: Sequence[np.ndarray]    # [E - E_min, E_max - E] in MeV
    dnde: np.ndarray                    # MeV-1 cm-2 s-1
    dnde_err: np.ndarray
    ul_energy: np.ndarray
    ul_energy_err: Sequence[np.ndarray]
    ul_e2dnde: np.ndarray               # MeV cm-2 s-1
    ul_yerr: np.ndarray                 # length of the UL arrows
    warning_mask: Optional[np.ndarray] = None  # bool, one per Fermi-LAT data point (less than 5 photons)
    n_VHE: int = 0
    n_VHE_UL: int = 0
    # EBL-corrected versions (filled by correct_for_EBL)
    dnde_EBL: Optional[np.ndarray] = None
    dnde_err_EBL: Optional[np.ndarray] = None
    ul_e2dnde_EBL: Optional[np.ndarray] = None
    ul_yerr_EBL: Optional[np.ndarray] = None

    @property
    def EBL_corrected(self):
        return self.dnde_EBL is not None

    @property
    def n_LAT(self):
        return len(self.energy) - self.n_VHE

    @property
    def n_LAT_UL(self):
        return len(self.ul_energy) - self.n_VHE_UL

    def correct_for_EBL(self, absorption, redshift):
        attenuation = EBL_attenuation(absorption, self.energy, redshift)
        self.dnde_EBL = self.dnde / attenuation
        self.dnde_err_EBL = self.dnde_err / attenuation
        attenuation_UL = EBL_attenuation(absorption, self.ul_energy, redshift)
        self.ul_e2dnde_EBL = self.ul_e2dnde / attenuation_UL
        self.ul_yerr_EBL = self.ul_yerr / attenuation_UL

    def fit_arrays(self):
        """(energy, dnde, dnde_err) used in the fit: EBL corrected, if available."""
        if self.EBL_corrected:
            return self.energy, self.dnde_EBL, self.dnde_err_EBL
        return self.energy, self.dnde, self.dnde_err

    def log_energy_grid(self, n=1000):
        """High-resolution grid in log10(E) covering all data points and upper limits."""
        xmin = np.log10(self.energy[0] - self.energy_err[0][0])
        xmax = np.log10(self.energy[-1] + self.energy_err[1][-1])
        if len(self.ul_energy) > 0:
            xmin = min(xmin, np.log10(self.ul_energy[0] - self.ul_energy_err[0][0]))
            xmax = max(xmax, np.log10(self.ul_energy[-1] + self.ul_energy_err[1][-1]))
        return np.linspace(xmin, xmax, n)

    @classmethod
    def from_sed_table(cls, sed, TSmin=TS_MIN):
        """
        Builds the data from a SED table (dict of arrays, as read by easyfermi_SED_in_sequence.read_SEDs,
        or the dict returned by gta.sed). Works with both e_ctr/e2dnde_ul95 (gta.sed dict) and
        e_ref/e2dnde_ul (fermipy *sed.fits, UL_CONF = 0.95).
        Data points with non-finite or non-positive flux/error are discarded (they cannot be fitted in log space).
        """
        e_min, e_max = _get_column(sed, "e_min"), _get_column(sed, "e_max")
        e_ctr = _get_column(sed, "e_ctr", "e_ref")
        if e_ctr is None:
            e_ctr = np.sqrt(e_min * e_max)
        e2dnde = _get_column(sed, "e2dnde")
        e2dnde_err = _get_column(sed, "e2dnde_err")
        e2dnde_ul = _get_column(sed, "e2dnde_ul95", "e2dnde_ul")
        ts = _get_column(sed, "ts")

        # Same fix used in compute_SED() for the empty TS column on macOS:
        scan = _get_column(sed, "dloglike_scan")
        if np.any(np.isnan(ts)) and scan is not None and scan.ndim == 2:
            ts = -2 * scan[:, 1]

        detected = ts > TSmin
        good = detected & np.isfinite(e2dnde) & np.isfinite(e2dnde_err) & (e2dnde > 0) & (e2dnde_err > 0)
        upper = (ts <= TSmin) & np.isfinite(e2dnde_ul)

        warning = _get_column(sed, "warning_few_photons")
        warning_mask = (warning[good] > 0) if warning is not None else None

        E = e_ctr[good]
        E_ul = e_ctr[upper]
        return cls(energy=E,
                   energy_err=[E - e_min[good], e_max[good] - E],
                   dnde=e2dnde[good] / E**2,
                   dnde_err=e2dnde_err[good] / E**2,
                   ul_energy=E_ul,
                   ul_energy_err=[E_ul - e_min[upper], e_max[upper] - E_ul],
                   ul_e2dnde=e2dnde_ul[upper],
                   ul_yerr=UL_YERR_FRACTION * e2dnde_ul[upper],
                   warning_mask=warning_mask)


# ----------------------------------------------------------------------------
# Likelihood
# ----------------------------------------------------------------------------
def _prior_box(model):
    lower = np.array([b[0] for b in model.bounds], dtype=float)
    upper = np.array([b[1] for b in model.bounds], dtype=float)
    return lower, upper


class SEDLikelihood:
    """
    Log-likelihood of a spectral model given a SED, vectorized over the walkers.

    Data points (TS > TSmin, Fermi-LAT and VHE, EBL corrected if requested): chi^2 in log10(dN/dE), with the errors
    propagated to log space (same as the previous easyfermi fit).

    Energy bins without detection (TS <= TSmin), only with method = "chi2_ul" (default):
        - if the bin has a best-fit value and error (e2dnde, e2dnde_err in the SED table), a gaussian term in linear
          flux, -0.5 * ((model - value) / error)^2. This is the parabolic approximation of the bin likelihood around its
          maximum: a model above the measured flux is penalized, a model below it is (almost) free. It uses only the
          best-fit normalization and its covariance error, not the likelihood scan (norm_scan/dloglike_scan), which is
          wrong in fermipy when the SED binning does not match binsperdec (and so are e2dnde_ul and the asymmetric errors);
        - otherwise (e.g. VHE upper limits, given only as upper limits), a censored gaussian term
          log Phi((UL - model) / sigma), with sigma = UL / 1.645, i.e. a 95% upper limit of a background-dominated
          measurement compatible with zero flux (Sawicki 2012, PASP 124, 1208).
    With method = "chi2", the upper limits are ignored (previous easyfermi behaviour).

    Parameters
    ----------
    model : SpectralModel
    log_e0 : float
        log10(pivot energy / MeV).
    data : SEDData
        Data points and upper limits (EBL corrected if data.correct_for_EBL was called). Its last data.n_VHE_UL upper
        limits are VHE upper limits.
    sed : dict-like or None
        Fermi-LAT SED table (gta.sed() output or *sed.fits table), used to get e2dnde and e2dnde_err of the bins without
        detection. If None, the Fermi-LAT upper limits of `data` enter as censored terms.
    absorption, redshift :
        EBL absorption model and redshift, used to correct the values of the bins without detection.
    """

    def __init__(self, model, log_e0, data, sed=None, absorption=None, redshift=0.0, method="chi2_ul", TSmin=TS_MIN):
        if method not in ("chi2_ul", "chi2"):
            raise ValueError("method must be 'chi2_ul' or 'chi2'.")
        self.model = model
        self.log_e0 = log_e0
        self.lower, self.upper = _prior_box(model)
        self.method = method

        # Data points: chi^2 in log space
        energy, dnde, dnde_err = data.fit_arrays()
        self.x_points = np.log10(energy)
        self.y_points = np.log10(dnde)                                  # log10(dN/dE)
        self.yerr_points = dnde_err / (dnde * np.log(10))

        # Bins without detection
        empty = np.zeros(0)
        self.x_gauss, self.value_gauss, self.err_gauss = empty, empty, empty   # gaussian terms (E^2 dN/dE)
        self.x_cens, self.ul_cens, self.sigma_cens = empty, empty, empty       # censored terms (E^2 dN/dE)
        if method == "chi2_ul":
            self._setup_upper_limits(data, sed, absorption, redshift, TSmin)

    def _setup_upper_limits(self, data, sed, absorption, redshift, TSmin):
        x_gauss, value_gauss, err_gauss, x_cens, ul_cens = [], [], [], [], []

        # Fermi-LAT bins without detection
        n_LAT_UL = data.n_LAT_UL
        lat_from_table = False
        if sed is not None:
            e_ctr = _get_column(sed, "e_ctr", "e_ref")
            if e_ctr is None:
                e_min, e_max = _get_column(sed, "e_min"), _get_column(sed, "e_max")
                e_ctr = np.sqrt(e_min * e_max) if e_min is not None and e_max is not None else None
            ts, value, error = _get_column(sed, "ts"), _get_column(sed, "e2dnde"), _get_column(sed, "e2dnde_err")
            ul = _get_column(sed, "e2dnde_ul95", "e2dnde_ul")
            if e_ctr is not None and ts is not None and value is not None and error is not None:
                lat_from_table = True
                no_detection = (ts <= TSmin) & np.isfinite(e_ctr)
                attenuation = np.ones(len(e_ctr))
                if absorption is not None and redshift > 0:
                    attenuation = EBL_attenuation(absorption, e_ctr, redshift)
                with_value = no_detection & np.isfinite(value) & np.isfinite(error) & (error > 0)
                x_gauss.append(np.log10(e_ctr[with_value]))
                value_gauss.append(value[with_value] / attenuation[with_value])
                err_gauss.append(error[with_value] / attenuation[with_value])
                if ul is not None:  # bins without value/error but with an upper limit
                    only_ul = no_detection & ~with_value & np.isfinite(ul) & (ul > 0)
                    x_cens.append(np.log10(e_ctr[only_ul]))
                    ul_cens.append(ul[only_ul] / attenuation[only_ul])
        if not lat_from_table and n_LAT_UL > 0:
            ul = data.ul_e2dnde_EBL if data.EBL_corrected else data.ul_e2dnde
            x_cens.append(np.log10(data.ul_energy[:n_LAT_UL]))
            ul_cens.append(ul[:n_LAT_UL])

        # VHE upper limits (only the upper limit is known)
        if data.n_VHE_UL > 0:
            ul = data.ul_e2dnde_EBL if data.EBL_corrected else data.ul_e2dnde
            x_cens.append(np.log10(data.ul_energy[n_LAT_UL:]))
            ul_cens.append(ul[n_LAT_UL:])

        if x_gauss:
            self.x_gauss, self.value_gauss, self.err_gauss = (np.concatenate(a) for a in (x_gauss, value_gauss, err_gauss))
        if x_cens:
            self.x_cens, self.ul_cens = np.concatenate(x_cens), np.concatenate(ul_cens)
            self.sigma_cens = self.ul_cens / 1.645

    @property
    def n_upper_limits(self):
        return len(self.x_gauss) + len(self.x_cens)

    # ---- evaluation (vectorized: thetas has shape (n_walkers, ndim)) -----------
    def _log_e2dnde(self, thetas, x):
        theta_columns = thetas.T[:, :, None]                                    # each parameter: (n_walkers, 1)
        return self.model.log_e2dnde(theta_columns, x, self.log_e0)           # (n_walkers, len(x))

    def log_likelihood(self, thetas):
        from scipy.special import log_ndtr

        thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
        loglike = np.zeros(len(thetas))
        with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
            if len(self.x_points) > 0:
                model = self._log_e2dnde(thetas, self.x_points) - 2 * self.x_points   # log10(dN/dE)
                loglike += -0.5 * np.sum(((self.y_points - model) / self.yerr_points)**2, axis=1)
            if len(self.x_gauss) > 0:
                model = 10 ** self._log_e2dnde(thetas, self.x_gauss)
                loglike += -0.5 * np.sum(((model - self.value_gauss) / self.err_gauss)**2, axis=1)
            if len(self.x_cens) > 0:
                model = 10 ** self._log_e2dnde(thetas, self.x_cens)
                loglike += np.sum(log_ndtr((self.ul_cens - model) / self.sigma_cens), axis=1)
        return np.where(np.isfinite(loglike), loglike, -np.inf)

    def log_probability(self, thetas):
        thetas = np.atleast_2d(np.asarray(thetas, dtype=float))
        inside = np.all((thetas > self.lower) & (thetas < self.upper), axis=1)  # uniform priors (also rejects NaN)
        result = np.full(len(thetas), -np.inf)
        if inside.any():
            result[inside] = self.log_likelihood(thetas[inside])
        return result

    def maximum_likelihood(self, n_prior_samples=N_PRIOR_SAMPLES, n_starts=5):
        """Maximum-likelihood point: random search over the prior box, refined with Nelder-Mead from the best points."""
        from scipy.optimize import minimize

        candidates = self.lower + (self.upper - self.lower) * np.random.rand(n_prior_samples, self.model.ndim)
        candidates = np.vstack([np.asarray(self.model.initial, dtype=float), candidates])
        values = self.log_probability(candidates)
        best_theta, best_value = candidates[np.argmax(values)], np.max(values)
        for start in candidates[np.argsort(values)[::-1][:n_starts]]:
            fit = minimize(lambda t: -self.log_probability(t)[0], start, method="Nelder-Mead",
                           options={"xatol": 1e-6, "fatol": 1e-6, "maxiter": 4000})
            if np.isfinite(fit.fun) and -fit.fun > best_value:
                best_theta, best_value = fit.x, -fit.fun
        return best_theta, best_value


def prepare_fit_data(model, energy, dnde, dnde_err):
    """Converts the data to the log10 space used in the chi^2 fit: (x, y, yerr)."""
    x = np.log10(energy)
    y = np.log10(dnde)
    yerr = dnde_err / (dnde * np.log(10))  # error propagation to log10 space (same for dnde and E^2 dnde)
    if model.fits_e2dnde:
        y = y + 2 * x
    return x, y, yerr


# ----------------------------------------------------------------------------
# MCMC
# ----------------------------------------------------------------------------
@dataclass
class MCMCResult:
    model: SpectralModel
    log_e0: float
    samples: np.ndarray      # flat chain after burn-in and thinning, shape (n_samples, ndim)
    log_prob: np.ndarray     # flat log-probability of the samples
    likelihood: SEDLikelihood
    theta_ml: np.ndarray     # maximum-likelihood point found by the optimizer
    diagnostics: dict        # likelihood method, walkers, steps, burn-in, thinning, tau, acceptance, convergence

    def __post_init__(self):
        # Highest-likelihood parameters (best of the optimizer and of the chain; the priors are uniform):
        candidates = np.vstack([self.theta_ml, self.samples[np.argmax(self.log_prob)]])
        loglike = self.likelihood.log_likelihood(candidates)
        self.theta_max = candidates[np.argmax(loglike)]
        self.max_loglike = float(np.max(loglike))
        self.AIC = 2 * self.model.ndim - 2 * self.max_loglike  # Akaike information criterion
        self.quantiles = np.quantile(self.samples, q=QUANTILES, axis=0).T  # shape (ndim, 3)

    # ---- summaries ---------------------------------------------------------
    def table_rows(self):
        """Rows (name, value, error_minus, error_plus) of the 'MCMC Parameters' table."""
        rows = [(name, q[1], q[1] - q[0], q[2] - q[1]) for name, q in zip(self.model.parameter_names, self.quantiles)]
        rows += [(name, value, 0.0, 0.0) for name, value in self.model.fixed_parameters]
        if self.model.uses_pivot:
            rows.append(("Ep=Emin (log scale)", self.log_e0, 0.0, 0.0))
        rows.append(("Akaike_IC", self.AIC, 0.0, 0.0))
        return rows

    def diagnostics_line(self):
        d = self.diagnostics
        return (f"MCMC diagnostics: likelihood = {d['likelihood']} ({d['n_upper_limits']} upper limits used), "
                f"walkers = {d['n_walkers']}, steps = {d['n_steps']} "
                f"(burn-in = {d['burn_in']}, thin = {d['thin']}), tau_max = {d['tau_max']:.1f}, "
                f"chain = {d['n_tau']:.1f} tau, acceptance = {d['acceptance']:.2f}, converged = {d['converged']}\n")

    def text_lines(self):
        """Lines written in Target_results.txt."""
        lines = [f"{name}: {q[1]} - {q[1] - q[0]} + {q[2] - q[1]}\n"
                 for name, q in zip(self.model.parameter_names, self.quantiles)]
        lines += [f"{name}: {value:g}\n" for name, value in self.model.fixed_parameters]
        lines.append(f"Akaike information criterion: {self.AIC}\n")
        lines.append(self.diagnostics_line())
        return lines

    def parameters_hdu(self):
        rows = self.table_rows()
        columns = [pyfits.Column(name="Parameter", array=np.array([r[0] for r in rows]), format="22A"),
                   pyfits.Column(name="Value", array=np.array([r[1] for r in rows], dtype=float), format="D"),
                   pyfits.Column(name="error_minus", array=np.array([r[2] for r in rows], dtype=float), format="D"),
                   pyfits.Column(name="error_plus", array=np.array([r[3] for r in rows], dtype=float), format="D")]
        hdu = pyfits.BinTableHDU.from_columns(columns, name="MCMC Parameters")
        d = self.diagnostics
        hdu.header["MODEL"] = (self.model.name, "Spectral model")
        hdu.header["LIKELIHD"] = (d["likelihood"], "chi2_ul: points + upper limits; chi2: points only")
        hdu.header["NUPLIM"] = (d["n_upper_limits"], "Bins without detection used in the fit")
        hdu.header["NWALKERS"] = (d["n_walkers"], "Number of walkers")
        hdu.header["NSTEPS"] = (d["n_steps"], "Steps per walker")
        hdu.header["NBURN"] = (d["burn_in"], "Burn-in steps discarded")
        hdu.header["NTHIN"] = (d["thin"], "Thinning of the saved posterior")
        hdu.header["TAUMAX"] = (round(d["tau_max"], 2), "Max. autocorrelation time (steps)")
        hdu.header["NTAU"] = (round(d["n_tau"], 2), "Chain length / TAUMAX")
        hdu.header["ACCEPT"] = (round(d["acceptance"], 3), "Mean acceptance fraction")
        hdu.header["CONVERGD"] = (bool(d["converged"]), "Convergence criterion reached")
        hdu.header["MAXLOGL"] = (self.max_loglike, "Maximum log-likelihood")
        return hdu

    def posterior_hdu(self):
        columns = [pyfits.Column(name=name, array=self.samples[:, i], format="D")
                   for i, name in enumerate(self.model.posterior_names)]
        return pyfits.BinTableHDU.from_columns(columns, name="MCMC Posterior dist.")

    # ---- model curves ------------------------------------------------------
    def e2dnde(self, x, theta=None):
        """E^2 dN/dE (MeV cm-2 s-1) of the model at x = log10(E/MeV). Default: highest-likelihood parameters."""
        theta = self.theta_max if theta is None else theta
        return 10 ** self.model.log_e2dnde(theta, x, self.log_e0)

    def random_samples(self, n=100):
        return self.samples[np.random.randint(len(self.samples), size=n)]


def run_MCMC(likelihood, nwalkers=N_WALKERS, chunk_steps=CHUNK_STEPS, n_tau_stop=N_TAU_STOP, tau_rtol=TAU_RTOL,
             burn_in_tau=BURN_IN_TAU, max_steps=MAX_STEPS, initial_spread=INITIAL_SPREAD, seed=None, verbose=True):
    """
    Adaptive MCMC (emcee, vectorized over the walkers).

    1) The maximum-likelihood point is found (random search over the prior + Nelder-Mead), and the walkers start
       in a small ball around it.
    2) The chain runs in chunks of `chunk_steps` steps. After each chunk the integrated autocorrelation time tau
       of every parameter is estimated. The run stops when the chain is longer than n_tau_stop * tau and tau changed
       by less than tau_rtol since the previous chunk (criterion of the emcee documentation), or at max_steps.
    3) The first burn_in_tau * tau steps are discarded and the chain is thinned by tau/2, so the saved posterior
       samples are nearly independent.

    Parameters
    ----------
    likelihood : SEDLikelihood
    seed : int or None
        Seed of numpy's global random generator (for reproducible results).

    Returns
    -------
    MCMCResult (convergence diagnostics in result.diagnostics).
    """
    if seed is not None:
        np.random.seed(seed)
    model = likelihood.model
    ndim = model.ndim
    nwalkers = max(int(nwalkers), 2 * ndim + 2)

    # 1) Starting point:
    theta_ml, _ = likelihood.maximum_likelihood()
    width = likelihood.upper - likelihood.lower
    p0 = theta_ml + initial_spread * width * np.random.randn(nwalkers, ndim)
    p0 = np.clip(p0, likelihood.lower + 1e-6 * width, likelihood.upper - 1e-6 * width)

    # 2) Adaptive run:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, likelihood.log_probability, vectorize=True)
    state, tau_old, converged = p0, None, False
    if verbose:
        print(f"Running MCMC ({model.name}, {likelihood.method} likelihood, {nwalkers} walkers)...")
    while sampler.iteration < max_steps:
        state = sampler.run_mcmc(state, min(chunk_steps, max_steps - sampler.iteration), progress=False)
        tau = sampler.get_autocorr_time(tol=0)
        if (tau_old is not None and np.all(np.isfinite(tau))
                and sampler.iteration > n_tau_stop * np.max(tau)
                and np.all(np.abs(tau_old - tau) / tau < tau_rtol)):
            converged = True
            break
        tau_old = tau

    # 3) Burn-in and thinning:
    tau_max = float(np.nanmax(tau)) if np.any(np.isfinite(tau)) else float(sampler.iteration)
    burn_in = int(min(burn_in_tau * tau_max, sampler.iteration // 2))
    thin = max(1, int(tau_max / 2))
    samples = sampler.get_chain(discard=burn_in, thin=thin, flat=True)
    log_prob = sampler.get_log_prob(discard=burn_in, thin=thin, flat=True)

    diagnostics = {"likelihood": likelihood.method, "n_upper_limits": likelihood.n_upper_limits, "n_walkers": nwalkers, "n_steps": int(sampler.iteration),
                   "burn_in": burn_in, "thin": thin, "tau_max": tau_max, "n_tau": sampler.iteration / tau_max,
                   "acceptance": float(np.mean(sampler.acceptance_fraction)), "converged": converged}
    if verbose:
        print(f"  {sampler.iteration} steps, tau_max = {tau_max:.1f} ({diagnostics['n_tau']:.0f} tau), "
              f"acceptance = {diagnostics['acceptance']:.2f}, {len(samples)} posterior samples.")
    if not converged:
        if diagnostics["n_tau"] >= MIN_N_TAU:
            print(f"- NOTE: the MCMC ({model.name}) stopped at {max_steps} steps with a chain of {diagnostics['n_tau']:.0f} "
                  f"autocorrelation times (< {n_tau_stop}): the quantiles are reliable, the estimate of tau less so.")
        else:
            print(f"- WARNING: the MCMC ({model.name}) did not converge: the chain has only {diagnostics['n_tau']:.1f} "
                  f"autocorrelation times after {max_steps} steps. The posterior is NOT reliable.")

    return MCMCResult(model=model, log_e0=likelihood.log_e0, samples=samples, log_prob=log_prob,
                      likelihood=likelihood, theta_ml=np.asarray(theta_ml, dtype=float), diagnostics=diagnostics)


# ----------------------------------------------------------------------------
# VHE data
# ----------------------------------------------------------------------------
def read_VHE_SED(filename, TSmin=TS_MIN):
    """
    Reads a VHE SED in the gammapy format (e_ref, e_min, e_max in TeV; e2dnde, e2dnde_err, e2dnde_ul in
    TeV cm-2 s-1; ts). Rows with NaN upper limit are discarded. Returns a dict (energies in MeV and fluxes in
    MeV cm-2 s-1), or None if the file cannot be read.
    """
    try:
        with pyfits.open(filename) as hdul:
            table = hdul[1].data
            keep = ~np.isnan(_as_float(table["e2dnde_ul"]))
            all_rows = {"energy": _as_float(table["e_ref"])[keep] * 1e6,      # TeV -> MeV
                        "energy_min": _as_float(table["e_min"])[keep] * 1e6,
                        "energy_max": _as_float(table["e_max"])[keep] * 1e6,
                        "e2dnde": _as_float(table["e2dnde"])[keep] * 1e6,     # TeV cm-2 s-1 -> MeV cm-2 s-1
                        "e2dnde_err": _as_float(table["e2dnde_err"])[keep] * 1e6,
                        "e2dnde_ul": _as_float(table["e2dnde_ul"])[keep] * 1e6,
                        "ts": _as_float(table["ts"])[keep]}
    except Exception as error:
        print(f"- WARNING: the VHE file could not be read ({type(error).__name__}: {error}). VHE data will not be used.")
        return None

    detected = all_rows["ts"] > TSmin
    upper = all_rows["ts"] <= TSmin
    E, E_ul = all_rows["energy"][detected], all_rows["energy"][upper]
    return {"all": all_rows,
            "energy": E,
            "energy_err": [E - all_rows["energy_min"][detected], all_rows["energy_max"][detected] - E],
            "e2dnde": all_rows["e2dnde"][detected],
            "e2dnde_err": all_rows["e2dnde_err"][detected],
            "ul_energy": E_ul,
            "ul_energy_err": [E_ul - all_rows["energy_min"][upper], all_rows["energy_max"][upper] - E_ul],
            "ul_e2dnde": all_rows["e2dnde_ul"][upper]}


def VHE_EBL_table(VHE, absorption, redshift):
    """HDU 'VHE data corrected for EBL' with all the VHE rows (NaN ULs removed)."""
    rows = VHE["all"]
    attenuation = EBL_attenuation(absorption, rows["energy"], redshift)
    columns = [pyfits.Column(name="energy", array=rows["energy"], format="D", unit="MeV"),
               pyfits.Column(name="energy_min", array=rows["energy_min"], format="D", unit="MeV"),
               pyfits.Column(name="energy_max", array=rows["energy_max"], format="D", unit="MeV"),
               pyfits.Column(name="e2dnde_VHE", array=rows["e2dnde"] / attenuation, format="D", unit="MeV cm-2 s-1"),
               pyfits.Column(name="e2dnde_VHE_err", array=rows["e2dnde_err"] / attenuation, format="D", unit="MeV cm-2 s-1"),
               pyfits.Column(name="e2dnde_VHE_UL95", array=rows["e2dnde_ul"] / attenuation, format="D", unit="MeV cm-2 s-1"),
               pyfits.Column(name="TS", array=rows["ts"], format="D")]
    return pyfits.BinTableHDU.from_columns(columns, name="VHE data corrected for EBL")


# ----------------------------------------------------------------------------
# FITS output
# ----------------------------------------------------------------------------
def EBL_corrected_SED_columns(e_ctr, e2dnde, e2dnde_err, e2dnde_ul95, absorption, redshift):
    """EBL-corrected columns to be added to the SED table (all energy bins)."""
    attenuation = EBL_attenuation(absorption, e_ctr, redshift)
    unit = "cm-2 MeV s-1"
    return [pyfits.Column(name="e2dnde_EBL_corrected", array=_as_float(e2dnde) / attenuation, format="D", unit=unit),
            pyfits.Column(name="e2dnde_err_EBL_corrected", array=_as_float(e2dnde_err) / attenuation, format="D", unit=unit),
            pyfits.Column(name="e2dnde_ul95_EBL_corrected", array=_as_float(e2dnde_ul95) / attenuation, format="D", unit=unit)]


def update_SED_fits(sed_file, new_sed_columns=None, new_hdus=()):
    """
    Adds columns to the SED table and appends HDUs to a *sed.fits file.
    Columns/HDUs with the same name already in the file are replaced, so the analysis can be re-run safely.
    """
    with pyfits.open(sed_file, memmap=False) as hdul:
        hdul.readall()

        if new_sed_columns:
            sed_index = next((i for i, h in enumerate(hdul) if h.name.upper() == "SED"), 1)
            new_names = {c.name.lower() for c in new_sed_columns}
            kept = [c for c in hdul[sed_index].columns if c.name.lower() not in new_names]
            table = pyfits.BinTableHDU.from_columns(pyfits.ColDefs(kept) + pyfits.ColDefs(list(new_sed_columns)))
            hdul[sed_index].data = table.data
            hdul[sed_index].name = "SED"

        for hdu in new_hdus:
            for i in reversed(range(1, len(hdul))):
                if hdul[i].name.upper() == hdu.name.upper():
                    del hdul[i]
            hdul.append(hdu)

        hdul.writeto(sed_file, overwrite=True)


# ----------------------------------------------------------------------------
# Plots
# ----------------------------------------------------------------------------
def corner_plot(result, filename):
    figure = corner.corner(result.samples, show_titles=True, labels=list(result.model.corner_labels),
                           plot_datapoints=True, quantiles=list(QUANTILES))
    figure.savefig(filename, bbox_inches="tight")
    plt.close(figure)


def _VHE_slice(n_total, n_VHE):
    return slice(n_total - n_VHE, n_total)


def _ylimits(data):
    """Y limits of the SED plot (same logic historically used by easyfermi)."""
    E2 = data.energy**2
    e2dnde, e2dnde_err = data.dnde * E2, data.dnde_err * E2
    if len(e2dnde) > 0:
        ymax = 2 * (e2dnde + e2dnde_err).max()
        ymin = 0.5 * (e2dnde - e2dnde_err).min()
        ymax = min(ymax, 4 * e2dnde.max())
        if data.EBL_corrected:
            ymax = max(ymax, 2 * np.max((data.dnde_EBL + data.dnde_err_EBL) * E2))
        ymin = max(ymin, e2dnde.min() / 5.0)
        if len(data.ul_e2dnde) > 0:
            yaux = (data.ul_e2dnde_EBL if data.EBL_corrected else data.ul_e2dnde).max()
            if yaux > ymax:
                ymax = 2 * yaux
            yaux = data.ul_e2dnde.min()
            if yaux < ymin:
                ymin = 0.5 * yaux
    else:
        ymax = 2 * data.ul_e2dnde.max()
        ymin = 0.5 * data.ul_e2dnde.min()
    return ymin, ymax


def plot_SED_MCMC(result, data, x, filename, model_energy=None, model_dnde=None, title="", n_curves=100,
                  figsize=(6, 5), dpi=250):
    """
    SED plot with the data, upper limits, fermipy model (if given), 100 random posterior models and the
    highest-likelihood MCMC model.
    """
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    ax.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax.yaxis.set_minor_locator(AutoMinorLocator(2))
    ax.tick_params(which="major", length=5, direction="in")
    ax.tick_params(which="minor", length=2.5, direction="in", bottom=True, top=True, left=True, right=True)
    ax.tick_params(bottom=True, top=True, left=True, right=True)

    E, E2 = data.energy, data.energy**2
    VHE = _VHE_slice(len(E), data.n_VHE)
    VHE_UL = _VHE_slice(len(data.ul_energy), data.n_VHE_UL)
    E_err_VHE = [data.energy_err[0][VHE], data.energy_err[1][VHE]]
    E_err_VHE_UL = [data.ul_energy_err[0][VHE_UL], data.ul_energy_err[1][VHE_UL]]
    point_style = dict(markeredgecolor="black", ecolor="black", zorder=130, fmt="o")

    # Fermi-LAT points with less than 5 photons:
    warning_energy = np.zeros(0)
    if data.warning_mask is not None and len(data.warning_mask) == data.n_LAT:
        warning_energy = E[:data.n_LAT][data.warning_mask]
    has_warnings = len(warning_energy) > 0

    alpha_original_data = 1
    if data.EBL_corrected:
        ax.errorbar(E, data.dnde_EBL * E2, xerr=data.energy_err, yerr=data.dnde_err_EBL * E2, color="C0",
                    label="Fermi-LAT EBL corrected", **point_style)
        ax.errorbar(data.ul_energy, data.ul_e2dnde_EBL, xerr=data.ul_energy_err, yerr=data.ul_yerr_EBL, uplims=True,
                    color="C0", **point_style)
        if has_warnings:
            ax.plot(warning_energy, data.dnde_EBL[:data.n_LAT][data.warning_mask] * warning_energy**2, "mo",
                    zorder=131, label="Less than 5 photons")
        if data.n_VHE > 0:
            ax.errorbar(E[VHE], data.dnde_EBL[VHE] * E2[VHE], xerr=E_err_VHE, yerr=data.dnde_err_EBL[VHE] * E2[VHE],
                        color="C1", label="VHE EBL corrected", **point_style)
        if data.n_VHE_UL > 0:
            ax.errorbar(data.ul_energy[VHE_UL], data.ul_e2dnde_EBL[VHE_UL], xerr=E_err_VHE_UL,
                        yerr=data.ul_yerr_EBL[VHE_UL], uplims=True, color="C1", **point_style)
        alpha_original_data = 0.4

    # Observed (non EBL-corrected) data:
    original_style = dict(alpha=alpha_original_data, zorder=120, fmt="o")
    ax.errorbar(E, data.dnde * E2, xerr=data.energy_err, yerr=data.dnde_err * E2, color="C0",
                label="Fermi-LAT", **original_style)
    if data.n_VHE > 0:
        ax.errorbar(E[VHE], data.dnde[VHE] * E2[VHE], xerr=E_err_VHE, yerr=data.dnde_err[VHE] * E2[VHE], color="C1",
                    label="VHE instrument", **original_style)
    if not data.EBL_corrected and has_warnings:
        ax.plot(warning_energy, data.dnde[:data.n_LAT][data.warning_mask] * warning_energy**2, "mo", zorder=131,
                label="Less than 5 photons")
    ax.errorbar(data.ul_energy, data.ul_e2dnde, xerr=data.ul_energy_err, yerr=data.ul_yerr, uplims=True,
                color="C0", **original_style)
    if data.n_VHE_UL > 0:
        ax.errorbar(data.ul_energy[VHE_UL], data.ul_e2dnde[VHE_UL], xerr=E_err_VHE_UL, yerr=data.ul_yerr[VHE_UL],
                    uplims=True, color="C1", **original_style)

    # Models:
    for theta in result.random_samples(n_curves):
        ax.plot(10**x, result.e2dnde(x, theta), color="r", zorder=0, alpha=0.1)  # posterior distribution
    if model_energy is not None and model_dnde is not None:
        model_energy, model_dnde = _as_float(model_energy), _as_float(model_dnde)
        ax.plot(model_energy, model_dnde * model_energy**2, color="gray", alpha=0.9, zorder=119, label="Fermipy fit")
    ax.plot(10**x, result.e2dnde(x), color="black", zorder=119, label="Highest Likelihood MCMC")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Energy [MeV]")
    ax.set_ylabel("E$^2dN/dE$ [MeV cm$^{-2}$ s$^{-1}$]")
    ax.set_title(title)
    ax.grid(which="both", linestyle=":")
    ax.set_ylim(*_ylimits(data))
    ax.set_xlim(0.8 * 10**x.min(), 1.2 * 10**x.max())
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig(filename, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------------------
# MCMC of the SEDs of every light-curve bin
# ----------------------------------------------------------------------------
def _column_name(parameter):
    return re.sub(r"\W+", "_", parameter).strip("_")


def _write_per_bin_summary(summary, model, filename):
    """ECSV table with one row per time bin (parameters, errors, AIC, status)."""
    from astropy.table import Table

    names = list(model.parameter_names)
    table = Table()
    for key in ("bin_index", "tmin_mjd", "tmax_mjd", "n_points", "n_upper_limits"):
        table[key] = [np.nan if s[key] is None else s[key] for s in summary]
    for name in names:
        column = _column_name(name)
        for k, suffix in enumerate(("", "_err_minus", "_err_plus")):
            table[column + suffix] = [s["parameters"][name][k] if s["parameters"] else np.nan for s in summary]
    table["AIC"] = [s["AIC"] for s in summary]
    table["MCMC_converged"] = [s.get("converged", False) for s in summary]
    table["MCMC_chain_length_tau"] = [s.get("n_tau", np.nan) for s in summary]
    table["status"] = [s["status"] for s in summary]
    table.meta["model"] = model.name
    table.write(filename, format="ascii.ecsv", overwrite=True)


def fit_SEDs_per_bin(SED_data, model_name, Emin, redshift=0.0, EBL_model=None, EBL_path=None, source_name="",
                     TSmin=TS_MIN, min_points=3, plot_format="png", make_corner_plots=False, write_summary=True,
                     likelihood_method="chi2_ul", mcmc_options=None, verbose=True):
    """
    Runs the MCMC on the SED of every time bin of a light curve.

    For each time bin:
        - the best-fit parameters ("MCMC Parameters") and posteriors ("MCMC Posterior dist.") are saved in the
          *sed.fits file (plus the EBL-corrected SED columns if redshift > 0);
        - a SED plot (Quickplot_SED_MCMC.<plot_format>) is saved in the time-bin directory.
    A summary table (MCMC_per_bin_summary_<model>.ecsv) is written in the light-curve directory.

    Parameters
    ----------
    SED_data : list of dict
        Output of easyfermi_SED_in_sequence.read_SEDs / SED_in_sequence (self.SED_per_bin_data in easyfermi).
    model_name : str
        One of MODELS.
    Emin : float
        Pivot energy (MeV) of the models with a pivot (the minimum energy of the analysis in easyfermi).
    redshift, EBL_model, EBL_path :
        If redshift > 0, the data are corrected for EBL absorption (EBL_model is the GUI label, see EBL_MODELS).
    TSmin : float
        Energy bins with TS > TSmin are data points; the others are upper limits.
    min_points : int
        Minimum number of data points needed to run the MCMC in a time bin.
    likelihood_method : {"chi2_ul", "chi2"}
        Likelihood used in the fit (see SEDLikelihood): data points + upper limits, or data points only.
    mcmc_options : dict or None
        Extra keyword arguments for run_MCMC (e.g. {"max_steps": 20000, "seed": 1}).

    Returns
    -------
    summary : list of dict (one per time bin) with bin_index, tmin_mjd, tmax_mjd, status, n_points,
        n_upper_limits, parameters ({name: (value, error_minus, error_plus)} or None), AIC, sed_file, plot_file.
    """
    model = get_model(model_name)
    log_e0 = np.log10(Emin)
    absorption = EBL_absorption_model(EBL_model, redshift, EBL_path) if redshift > 0.0 else None

    summary = []
    for n, entry in enumerate(SED_data):
        bin_label = os.path.basename(os.path.normpath(entry["bin_directory"]))
        info = {"bin_index": entry.get("bin_index"), "tmin_mjd": entry.get("tmin_mjd"),
                "tmax_mjd": entry.get("tmax_mjd"), "sed_file": entry.get("sed_file"), "plot_file": None,
                "status": "", "n_points": 0, "n_upper_limits": 0, "parameters": None, "AIC": np.nan,
                "converged": False, "n_tau": np.nan}
        try:
            data = SEDData.from_sed_table(entry["sed"], TSmin)
            info["n_points"], info["n_upper_limits"] = len(data.energy), len(data.ul_energy)
            if len(data.energy) < min_points:
                info["status"] = f"skipped (only {len(data.energy)} data points with TS > {TSmin:g})"
                continue

            if absorption is not None:
                data.correct_for_EBL(absorption, redshift)
            if verbose:
                print(f"MCMC of {bin_label} ({n + 1}/{len(SED_data)}):")
            likelihood = SEDLikelihood(model, log_e0, data, sed=entry["sed"], absorption=absorption,
                                       redshift=redshift, method=likelihood_method, TSmin=TSmin)
            result = run_MCMC(likelihood, verbose=verbose, **(mcmc_options or {}))

            # Saving the results in the *sed.fits file of the time bin:
            new_columns = None
            if absorption is not None:
                sed = entry["sed"]
                e_ctr = _get_column(sed, "e_ctr", "e_ref")
                new_columns = EBL_corrected_SED_columns(e_ctr, _get_column(sed, "e2dnde"),
                                                        _get_column(sed, "e2dnde_err"),
                                                        _get_column(sed, "e2dnde_ul95", "e2dnde_ul"),
                                                        absorption, redshift)
            update_SED_fits(entry["sed_file"], new_columns, [result.parameters_hdu(), result.posterior_hdu()])

            # Plots:
            title = f"{source_name} - {model.name}"
            if entry.get("tmin_mjd") is not None:
                title = f"{source_name} - MJD {entry['tmin_mjd']:.2f}-{entry['tmax_mjd']:.2f} - {model.name}"
            model_flux = entry.get("model_flux") or {}
            info["plot_file"] = os.path.join(entry["bin_directory"], f"Quickplot_SED_MCMC.{plot_format}")
            plot_SED_MCMC(result, data, data.log_energy_grid(), info["plot_file"],
                          model_energy=_get_column(model_flux, "energy", "energies") if model_flux else None,
                          model_dnde=_get_column(model_flux, "dnde") if model_flux else None, title=title)
            if make_corner_plots:
                corner_plot(result, os.path.join(entry["bin_directory"], "Quickplot_MCMC_SED_pars.png"))

            info["parameters"] = {name: (q[1], q[1] - q[0], q[2] - q[1])
                                  for name, q in zip(model.parameter_names, result.quantiles)}
            info["AIC"] = float(result.AIC)
            info["converged"] = bool(result.diagnostics["converged"])
            info["n_tau"] = float(result.diagnostics["n_tau"])
            info["status"] = "done"

        except Exception as error:
            info["status"] = f"failed ({type(error).__name__}: {error})"
        finally:
            summary.append(info)
            if verbose:
                print(f"  {bin_label}: {info['status']}")

    n_done = sum(s["status"] == "done" for s in summary)
    print(f"MCMC per time bin: {n_done} done, {len(summary) - n_done} skipped/failed.")

    if write_summary and len(SED_data) > 0:
        LC_directory = os.path.dirname(os.path.normpath(SED_data[0]["bin_directory"]))
        filename = os.path.join(LC_directory, f"MCMC_per_bin_summary_{model.name}.ecsv")
        try:
            _write_per_bin_summary(summary, model, filename)
            print(f"Summary table saved in {filename}")
        except Exception as error:
            print(f"- WARNING: the summary table could not be written ({type(error).__name__}: {error}).")

    return summary
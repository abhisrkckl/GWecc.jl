"""Python interface for GWecc.jl to be used with ENTERPRISE.
Provides the PTA signal due to an eccentric supermassive binary
as an ENTERPRISE Signal object."""

__version__ = "0.1.1"

import json

import enterprise
import numpy as np
from enterprise.signals.deterministic_signals import Deterministic
from enterprise.signals.parameter import (
    Parameter,
    TruncNormalPrior,
    TruncNormalSampler,
    Uniform,
)
from enterprise.signals.signal_base import function as enterprise_function
from juliacall import Main as jl
from scipy.special import erf

jl.seval("using GWecc")


@enterprise_function
def eccentric_pta_signal_1psr(
    toas,
    sigma,
    rho,
    log10_M,
    eta,
    log10_F,
    e0,
    l0,
    tref,
    log10_A,
    deltap,
    psrTerm=False,
    spline=False,
):
    """Compute the eccentric SMBHB PTA signal for the single-pulsar case
    using the reduced parametrization. This is a thin wrapper around
    `GWecc.eccentric_pta_signal_1psr`.

    Ref: Susobhanan 2023

    Parameters
    ----------
    toas : array-like
        Collection of TOAs (s)
    sigma : enterprise.signals.parameter.Parameter
        Projection angle 1 (rad)
    rho : enterprise.signals.parameter.Parameter
        Projection angle 2 (rad)
    log10_M : enterprise.signals.parameter.Parameter
        Log10 total mass (Msun)
    eta : enterprise.signals.parameter.Parameter
        Symmetric mass ratio
    log10_F : enterprise.signals.parameter.Parameter
        Log10 initial GW frequency (Hz)
    e0 : enterprise.signals.parameter.Parameter
        Initial eccentricity
    l0 : enterprise.signals.parameter.Parameter
        Initial mean anomaly (rad)
    tref : float
        Fiducial time (s)
    log10_A : enterprise.signals.parameter.Parameter
        Log10 effective PTA signal amplitude (s)
    deltap : enterprise.signals.parameter.Parameter
        Pulsar term delay (yr)
    psrTerm : bool
        Whether to include pulsar term (default is False)
    spline : bool
        Whether to use spline-based fast computation (default is False)

    Returns
    -------
    Rs : array-like
        The PTA signal (s) evaluated at TOAs.

    Notes
    -----
    1. This thin wrapper is required because ENTERPRISE relies on reflection
       of Python functions, which does not work properly with juliacall.

    2. This is an `enterprise_function`, and can be partially evaluated.

    3. deltap will be ignored if psrTerm is False.
    """
    return jl.eccentric_pta_signal_1psr(
        toas,
        float(sigma),
        float(rho),
        float(log10_M),
        float(eta),
        float(log10_F),
        float(e0),
        float(l0),
        float(tref),
        float(log10_A),
        float(deltap),
        psrTerm,
        spline,
    )


@enterprise_function
def eccentric_pta_signal(
    toas,
    theta,
    phi,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    log10_M,
    eta,
    log10_F,
    e0,
    gamma0,
    gammap,
    l0,
    lp,
    tref,
    log10_A,
    psrdist,
    psrTerm=False,
    spline=False,
):
    """Compute the eccentric SMBHB PTA signal. This is a thin wrapper around
    `GWecc.eccentric_pta_signal`.

    Ref: Susobhanan 2023

    Parameters
    ----------
    toas : array-like
        Collection of TOAs (s)
    theta : float
        Pulsar zenith angle (rad)
    phi : float
        Pulsar right ascension (rad)
    cos_gwtheta : float
        Cos zenith angle of the GW source (rad)
    gwphi : float
        Right ascension of the GW source (rad)
    psi : float
        Polarization angle (rad)
    cos_inc : float
        Cos inclination
    log10_M : float
        Log10 total mass (Msun)
    eta : float
        Symmetric mass ratio
    log10_F : float
        Log10 initial GW frequency (Hz)
    e0 : float
        Initial eccentricity
    gamma0 : float
        Initial periastron angle for the Earth term (rad)
    gammap : float
        Initial periastron angle for the Pulsar term (rad)
    l0 : float
        Initial mean anomaly for the Earth term (rad)
    lp : float
        Initial mean anomaly for the Pulsar term (rad)
    tref : float
        Fiducial time (s)
    log10_A : float
        Log10 PTA signal amplitude (s)
    psrdist : float
        Pulsar distance (kpc)
    psrTerm : bool
        Whether to include the pulsar term (default is False)
    spline : bool
        Whether to use the spline method for fast evaluation (default is False)

    Returns
    -------
    Rs : array-like
        The PTA signal (s) evaluated at TOAs.

    Notes
    -----
    1. This thin wrapper is required because ENTERPRISE relies on reflection
       of Python functions, which does not work properly with juliacall.

    2. This is an `enterprise_function`, and can be partially evaluated.

    3. lp, gammap, and delta_pdist will be ignored if psrTerm is False.
    """

    return jl.eccentric_pta_signal(
        toas,
        float(theta),
        float(phi),
        float(psrdist),
        float(cos_gwtheta),
        float(gwphi),
        float(psi),
        float(cos_inc),
        float(log10_M),
        float(eta),
        float(log10_F),
        float(e0),
        float(gamma0),
        float(gammap),
        float(l0),
        float(lp),
        float(tref),
        float(log10_A),
        psrTerm,
        spline,
    )


@enterprise_function
def eccentric_pta_signal_target(
    toas,
    theta,
    phi,
    cos_gwtheta,
    gwphi,
    psi,
    cos_inc,
    eta,
    log10_F,
    e0,
    gamma0,
    gammap,
    l0,
    lp,
    tref,
    log10_A,
    gwdist,
    psrdist,
    psrTerm=False,
    spline=False,
):
    """Compute the eccentric SMBHB PTA signal. This is a thin wrapper around
    `GWecc.eccentric_pta_signal_target`.

    Ref: Susobhanan 2023

    Parameters
    ----------
    toas : array-like
        Collection of TOAs (s)
    theta : float
        Pulsar zenith angle (rad)
    phi : float
        Pulsar right ascension (rad)
    cos_gwtheta : float
        Cos zenith angle of the GW source (rad)
    gwphi : float
        Right ascension of the GW source (rad)
    psi : float
        Polarization angle (rad)
    cos_inc : float
        Cos inclination
    eta : float
        Symmetric mass ratio
    log10_F : float
        Log10 initial GW frequency (Hz)
    e0 : float
        Initial eccentricity
    gamma0 : float
        Initial periastron angle for the Earth term (rad)
    gammap : float
        Initial periastron angle for the Pulsar term (rad)
    l0 : float
        Initial mean anomaly for the Earth term (rad)
    lp : float
        Initial mean anomaly for the Pulsar term (rad)
    tref : float
        Fiducial time (s)
    log10_A : float
        Log10 PTA signal amplitude (s)
    gwdist : float
        Luminosity distance to the GW source (Mpc)
    psrdist : float
        Pulsar distance (kpc)
    psrTerm : bool
        Whether to include the pulsar term (default is False)
    spline : bool
        Whether to use the spline method for fast evaluation (default is False)

    Returns
    -------
    Rs : array-like
        The PTA signal (s) evaluated at TOAs.

    Notes
    -----
    1. This thin wrapper is required because ENTERPRISE relies on reflection
       of Python functions, which does not work properly with juliacall.

    2. This is an `enterprise_function`, and can be partially evaluated.

    3. lp, gammap, and delta_pdist will be ignored if psrTerm is False.

    4. The total mass of the GW source is computed using the amplitude and the
       luminosity distance.
    """

    return jl.eccentric_pta_signal_target(
        toas,
        float(theta),
        float(phi),
        float(psrdist),
        float(cos_gwtheta),
        float(gwphi),
        float(psi),
        float(cos_inc),
        float(eta),
        float(log10_F),
        float(e0),
        float(gamma0),
        float(gammap),
        float(l0),
        float(lp),
        float(tref),
        float(log10_A),
        float(gwdist),
        psrTerm,
        spline,
    )


def gwecc_1psr_block(
    tref,
    sigma=Uniform(0, np.pi)("gwecc_sigma"),
    rho=Uniform(-np.pi, np.pi)("gwecc_rho"),
    log10_M=Uniform(6, 9)("gwecc_log10_M"),
    eta=Uniform(0, 0.25)("gwecc_eta"),
    log10_F=Uniform(-9, -7)("gwecc_log10_F"),
    e0=Uniform(0.01, 0.8)("gwecc_e0"),
    l0=Uniform(0, 2 * np.pi)("gwecc_l0"),
    log10_A=Uniform(-11, -7)("gwecc_log10_A"),
    deltap=Uniform(0, 2000),
    psrTerm=False,
    spline=False,
    name="gwecc",
):
    """Deterministic eccentric-orbit continuous GW model for a single pulsar
    using a reduced parametrization to avoid degeneracies. This should not be used
    while analyzing more than one pulsar.

    Ref: Susobhanan 2023

    Parameters
    ----------
    tref : float
        Fiducial time (s)
    sigma : enterprise.signals.parameter.Parameter
        Projection angle 1 (rad)
    rho : enterprise.signals.parameter.Parameter
        Projection angle 2 (rad)
    log10_M : enterprise.signals.parameter.Parameter
        Log10 total mass (Msun)
    eta : enterprise.signals.parameter.Parameter
        Symmetric mass ratio
    log10_F : enterprise.signals.parameter.Parameter
        Log10 initial GW frequency (Hz)
    e0 : enterprise.signals.parameter.Parameter
        Initial eccentricity
    l0 : enterprise.signals.parameter.Parameter
        Initial mean anomaly (rad)
    log10_A : enterprise.signals.parameter.Parameter
        Log10 effective PTA signal amplitude (s)
    deltap : enterprise.signals.parameter.Parameter
        Pulsar term delay (yr)
    psrTerm : bool
        Whether to include the pulsar term (default is False)
    spline : bool
        Whether to use spline-based fast computation (default is False)
    name : str
        Name of the signal object (default is "gwecc")

    Returns
    -------
    wf : enterprise.signals.deterministic_signals.Deterministic
        ENTERPRISE deterministic signal

    Notes
    -----
    1. deltap will be ignored if psrTerm is False.
    """

    deltap = deltap if psrTerm else 0.0

    return Deterministic(
        eccentric_pta_signal_1psr(
            sigma=sigma,
            rho=rho,
            log10_M=log10_M,
            eta=eta,
            log10_F=log10_F,
            e0=e0,
            l0=l0,
            tref=tref,
            log10_A=log10_A,
            deltap=deltap,
            psrTerm=psrTerm,
            spline=spline,
        ),
        name=name,
    )


def gwecc_block(
    tref,
    cos_gwtheta=Uniform(-1, 1)("gwecc_cos_gwtheta"),
    gwphi=Uniform(0, 2 * np.pi)("gwecc_gwphi"),
    psi=Uniform(0, np.pi)("gwecc_psi"),
    cos_inc=Uniform(-1, 1)("gwecc_cos_inc"),
    log10_M=Uniform(6, 9)("gwecc_log10_M"),
    eta=Uniform(0, 0.25)("gwecc_eta"),
    log10_F=Uniform(-9, -7)("gwecc_log10_F"),
    e0=Uniform(0.01, 0.8)("gwecc_e0"),
    gamma0=Uniform(0, np.pi)("gwecc_gamma0"),
    gammap=Uniform(0, np.pi),
    l0=Uniform(0, 2 * np.pi)("gwecc_l0"),
    lp=Uniform(0, 2 * np.pi),
    log10_A=Uniform(-11, -7)("gwecc_log10_A"),
    psrdist=None,
    psrTerm=False,
    tie_psrTerm=False,
    spline=False,
    name="gwecc",
):
    """Deterministic eccentric-orbit continuous GW model.

    Ref: Susobhanan 2023

    Parameters
    ----------
    tref : float
        Fiducial time (s)
    cos_gwtheta : enterprise.signals.parameter.Parameter
        Cos zenith angle of the GW source (rad)
    gwphi : enterprise.signals.parameter.Parameter
        Right ascension of the GW source (rad)
    psi : enterprise.signals.parameter.Parameter
        Polarization angle
    cos_inc : enterprise.signals.parameter.Parameter
        Cos inclination
    log10_M : enterprise.signals.parameter.Parameter
        Log10 total mass (Msun)
    eta : enterprise.signals.parameter.Parameter
        Symmetric mass ratio
    log10_F : enterprise.signals.parameter.Parameter
        Log10 initial GW frequency (Hz)
    e0 : enterprise.signals.parameter.Parameter
        Initial eccentricity
    gamma0 : enterprise.signals.parameter.Parameter
        Initial periastron angle of the Earth term (rad)
    gammap : enterprise.signals.parameter.Parameter
        Initial periastron angle of the Pulsar term (rad)
    l0 : enterprise.signals.parameter.Parameter
        Initial mean anomaly of the Earth term (rad)
    lp : enterprise.signals.parameter.Parameter
        Initial mean anomaly of the Pulsar term (rad)
    log10_A : enterprise.signals.parameter.Parameter
        Log10 effective PTA signal amplitude (s)
    psrdist : enterprise.signals.parameter.Parameter
        Pulsar distance (kpc)
    psrTerm : bool
        Whether to include the pulsar term (default is False)
    tie_psrTerm : bool
        Whether to tie the Pulsar term phase (lp, gammap) with
        the Earth term phase (l0, gamma0)
    spline : bool
        Whether to use spline-based fast computation (default is False)
    name : str
        Name of the signal object (default is "gwecc")

    Returns
    -------
    wf : enterprise.signals.deterministic_signals.Deterministic
        ENTERPRISE deterministic signal

    Notes
    -----
    1. lp, gammap, and delta_pdist will be ignored if psrTerm is False.

    2. lp and gammap will be ignored if tie_psrTerm is True.

    3. delta_pdist should have a truncated normal distribution to avoid
       negative values. Negative values will cause exceptions in the orbital
       dynamics code.
    """

    if psrTerm:
        assert (
            psrdist is not None and psrdist != 0
        ), "psrdist should be provided when psrTerm is True."
        if tie_psrTerm:
            # Tie pulsar term phase to the earth term phase.
            # This ignores the explicitly given gammap and lp
            gammap, lp = (gamma0, l0)
    else:
        # These are not used.
        gammap, lp = (0.0, 0.0)
        psrdist = 1.0

    return Deterministic(
        eccentric_pta_signal(
            cos_gwtheta=cos_gwtheta,
            gwphi=gwphi,
            psi=psi,
            cos_inc=cos_inc,
            log10_M=log10_M,
            eta=eta,
            log10_F=log10_F,
            e0=e0,
            gamma0=gamma0,
            gammap=gammap,
            l0=l0,
            lp=lp,
            tref=tref,
            log10_A=log10_A,
            psrdist=psrdist,
            psrTerm=psrTerm,
            spline=spline,
        ),
        name=name,
    )


def gwecc_target_block(
    tref,
    cos_gwtheta,
    gwphi,
    gwdist,
    psi=Uniform(0, np.pi)("gwecc_psi"),
    cos_inc=Uniform(-1, 1)("gwecc_cos_inc"),
    eta=Uniform(0, 0.25)("gwecc_eta"),
    log10_F=Uniform(-9, -7)("gwecc_log10_F"),
    e0=Uniform(0.01, 0.8)("gwecc_e0"),
    gamma0=Uniform(0, np.pi)("gwecc_gamma0"),
    gammap=Uniform(0, np.pi),
    l0=Uniform(0, 2 * np.pi)("gwecc_l0"),
    lp=Uniform(0, 2 * np.pi),
    log10_A=Uniform(-11, -7)("gwecc_log10_A"),
    psrdist=None,
    psrTerm=False,
    tie_psrTerm=False,
    spline=False,
    name="gwecc",
):
    """Deterministic eccentric-orbit continuous GW model for targeted
    search (known source location and distance).

    Ref: Susobhanan 2023

    Parameters
    ----------
    tref : float
        Fiducial time (s)
    cos_gwtheta : float
        Cos zenith angle of the GW source (rad)
    gwphi : float
        Right ascension of the GW source (rad)
    gwdist : float
        Luminosity distance to the GW source (Mpc)
    psi : enterprise.signals.parameter.Parameter
        Polarization angle (rad)
    cos_inc : enterprise.signals.parameter.Parameter
        Cos inclination
    eta : enterprise.signals.parameter.Parameter
        Symmetric mass ratio
    log10_F : enterprise.signals.parameter.Parameter
        Log10 initial GW frequency (Hz)
    e0 : enterprise.signals.parameter.Parameter
        Initial eccentricity
    gamma0 : enterprise.signals.parameter.Parameter
        Initial periastron angle of the Earth term (rad)
    gammap : enterprise.signals.parameter.Parameter
        Initial periastron angle of the Pulsar term (rad)
    l0 : enterprise.signals.parameter.Parameter
        Initial mean anomaly of the Earth term (rad)
    lp : enterprise.signals.parameter.Parameter
        Initial mean anomaly of the Pulsar term (rad)
    log10_A : enterprise.signals.parameter.Parameter
        Log10 effective PTA signal amplitude (s)
    psrdist : enterprise.signals.parameter.Parameter
        Pulsar distance (kpc)
    psrTerm : bool
        Whether to include the pulsar term (default is False)
    tie_psrTerm : bool
        Whether to tie the Pulsar term phase (lp, gammap) with
        the Earth term phase (l0, gamma0)
    spline : bool
        Whether to use spline-based fast computation (default is False)
    name : str
        Name of the signal object (default is "gwecc")

    Returns
    -------
    wf : enterprise.signals.deterministic_signals.Deterministic
        ENTERPRISE deterministic signal

    Notes
    -----
    1. lp, gammap, and delta_pdist will be ignored if psrTerm is False.

    2. lp and gammap will be ignored if tie_psrTerm is True.

    3. delta_pdist should have a truncated normal distribution to avoid
       negative values. Negative values will cause exceptions in the orbital
       dynamics code.

    4. The total mass of the GW source is computed using the amplitude and the
       luminosity distance.
    """

    if psrTerm:
        assert (
            psrdist is not None and psrdist != 0
        ), "psrdist should be provided when psrTerm is True."
        if tie_psrTerm:
            # Tie pulsar term phase to the earth term phase.
            # This ignores the explicitly given gammap and lp
            gammap, lp = (gamma0, l0)
    else:
        # These are not used.
        gammap, lp = (0.0, 0.0)
        psrdist = 1.0

    return Deterministic(
        eccentric_pta_signal_target(
            cos_gwtheta=cos_gwtheta,
            gwphi=gwphi,
            psi=psi,
            cos_inc=cos_inc,
            eta=eta,
            log10_F=log10_F,
            e0=e0,
            gamma0=gamma0,
            gammap=gammap,
            l0=l0,
            lp=lp,
            tref=tref,
            log10_A=log10_A,
            psrdist=psrdist,
            gwdist=gwdist,
            psrTerm=psrTerm,
            spline=spline,
        ),
        name=name,
    )


def gwecc_prior(pta, tref, tmax, name="gwecc"):
    def gwecc_target_prior_fn(params):
        param_map = pta.map_params(params)
        if jl.validate_params(
            param_map[f"{name}_log10_M"],
            param_map[f"{name}_eta"],
            param_map[f"{name}_log10_F"],
            param_map[f"{name}_e0"],
            tref,
            tmax,
        ):
            return pta.get_lnprior(param_map)
        else:
            return -np.inf

    return gwecc_target_prior_fn


def gwecc_target_prior(pta, gwdist, tref, tmax, log10_F=None, name="gwecc"):
    def gwecc_target_prior_fn(params):
        param_map = pta.map_params(params)

        log10_F = (
            param_map[f"{name}_log10_F"] if f"{name}_log10_F" in param_map else log10_F
        )
        assert (
            log10_F is not None
        ), "log10_F should either be given while calling gwecc_target_prior or be a model parameter."

        if jl.validate_params_target(
            param_map[f"{name}_log10_A"],
            param_map[f"{name}_eta"],
            log10_F,
            param_map[f"{name}_e0"],
            gwdist,
            tref,
            tmax,
        ):
            return pta.get_lnprior(param_map)
        else:
            return -np.inf

    return gwecc_target_prior_fn


def PsrDistPrior(psrdist_info: dict, dmdist_broaden_factor=2):
    """Truncated normal distribution for pulsar distances (psrdist).
    This is a truncated Normal distribution with mean and variance
    taken from the distance measurements, and a lower cutoff at 0, so
    that the pulsar distance doesn't go negative.

    The pulsar name is obtained by parsing the parameter name.
    The pulsar distance info is obtained from the following sources
    in that order:
        1. The `psrdist_info` dictionary
        2. The `pulsar_distances.json` file available in the ENTERPRISE
           distribution
        3. The default value used in ENTERPRISE (1.0 +/- 0.2)
    """

    class PsrDistPrior(Parameter):
        _pulsar_distance_info = psrdist_info
        _dmdist_broaden_factor = dmdist_broaden_factor
        _typename = "PsrDistPrior"
        _size = None

        def __init__(self, name):
            super().__init__(name)

            self.psrname = name.split("_")[0]

            if self.psrname in self._pulsar_distance_info:
                (
                    self.pdist,
                    self.pdist_sigma,
                    self.pdist_type,
                ) = self._pulsar_distance_info[self.psrname]
                assert self.pdist_type in ["PX", "DM"]

                if self.pdist_type == "DM":
                    # Broader distribution to account for systematic uncertainties
                    # in DM distance measurement
                    self.pdist_sigma *= self._dmdist_broaden_factor
            else:
                # Try with ENTERPRISE
                pdist_file_default = (
                    f"{enterprise.__path__[0]}/datafiles/pulsar_distances.json"
                )
                with open(pdist_file_default, "r") as f:
                    pdist_info_default = json.load(f)

                if self.psrname in pdist_info_default:
                    self.pdist, self.pdist_sigma = pdist_info_default[self.psrname]
                else:
                    # ENTERPRISE default value
                    self.pdist, self.pdist_sigma = 1, 0.2

            # Precompute the norm because erf is slow
            self.pnorm = (
                np.sqrt(2 / np.pi)
                / self.pdist_sigma
                / (1 + erf(self.pdist / np.sqrt(2) / self.pdist_sigma))
            )

        def _prior(self, name):
            """Hack to make sure Parameter.__init__ doesn't throw an error"""
            pass

        def prior(self, x):
            return TruncNormalPrior(
                x, self.pdist, self.pdist_sigma, 0, np.inf, norm=self.pnorm
            )

        def sample(self):
            return TruncNormalSampler(
                self.pdist, self.pdist_sigma, 0, np.inf, norm=self.pnorm
            )

        def __repr__(self):
            return self._typename

        @property
        def params(self):
            """Hack to make sure ENTERPRISE doesn't throw an error"""
            return [self]

    return PsrDistPrior

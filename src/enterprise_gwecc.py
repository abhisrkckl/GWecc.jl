"""Python interface for GWecc.jl to be used with ENTERPRISE.
Provides the PTA signal due to an eccentric supermassive binary
as an ENTERPRISE Signal object."""

__version__ = "0.1.0"

import numpy as np
from enterprise.signals.deterministic_signals import Deterministic
from enterprise.signals.parameter import Uniform
from enterprise.signals.signal_base import function as enterprise_function
from juliacall import Main as jl

jl.seval("using GWecc")

# This thin wrapper function is required because ENTERPRISE relies on reflection
# of Python functions, which does not work properly with juliacall.
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


# This thin wrapper function is required because ENTERPRISE relies on reflection
# of Python functions, which does not work properly with juliacall.
@enterprise_function
def eccentric_pta_signal(
    toas,
    theta,
    phi,
    pdist,
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
    psrTerm=False,
    spline=False,
):
    return jl.eccentric_pta_signal(
        toas,
        float(theta),
        float(phi),
        float(pdist[0] if isinstance(pdist, tuple) else pdist),
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


def gwecc_1psr_block(
    tref,
    sigma=Uniform(0, np.pi)("gwecc_sigma"),
    rho=Uniform(-1, 1)("gwecc_rho"),
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
    using a reduced parametrization to avoid degeneracies. Should not be used
    while analyzing more than one pulsar.

    Parameters
    ----------
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
        Log10 PTA signal amplitude (s)
    deltap : enterprise.signals.parameter.Parameter
        Pulsar term delay (yr)
    psrTerm : bool
        Whether to include pulsar term
    spline : bool
        Whether to use spline-based fast computation
    name : str
        Name of the signal object

    Returns
    -------
        ENTERPRISE deterministic signal (`enterprise.signals.deterministic_signals.Deterministic`)
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
    psrTerm=False,
    spline=False,
    name="gwecc",
):
    """Deterministic eccentric-orbit continuous GW model."""

    gammap, lp = (gammap, lp) if psrTerm else (0.0, 0.0)

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
            psrTerm=psrTerm,
            spline=spline,
        ),
        name=name,
    )

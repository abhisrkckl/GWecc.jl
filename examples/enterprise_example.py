"""This example demonstrates how to create an ENTERPRISE Signal object from GWecc."""

from enterprise.signals.signal_base import function as enterprise_function
from enterprise.signals.deterministic_signals import Deterministic
from enterprise.signals.parameter import Uniform

import numpy as np

from juliacall import Main as jl

jl.seval("using GWecc")

# This thin wrapper function is required because ENTERPRISE relies on reflection
# of Python functions, which does not work properly with juliacall.
@enterprise_function
def eccentric_pta_signal_planck18_1psr(
    toas,
    pdist,
    alpha,
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
    log10_zc,
    psrTerm=False,
):
    return jl.eccentric_pta_signal_planck18_1psr(
        toas,
        pdist,
        alpha,
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
        log10_zc,
        psrTerm,
    )


# This thin wrapper function is required because ENTERPRISE relies on reflection
# of Python functions, which does not work properly with juliacall.
@enterprise_function
def eccentric_pta_signal_planck18(
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
    log10_zc,
    psrTerm=False,
):
    return jl.eccentric_pta_signal_planck18(
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
        log10_zc,
        psrTerm,
    )


def gwecc_1psr_block(
    tref,
    alpha=Uniform(0, 1)("gwecc_alpha"),
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
    log10_zc=Uniform(-4, -3)("gwecc_log10_zc"),
):
    return Deterministic(
        eccentric_pta_signal_planck18_1psr(
            alpha=alpha,
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
            log10_zc=log10_zc,
        ),
        name="gwecc",
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
    log10_zc=Uniform(-4, -3)("gwecc_log10_zc"),
):
    return Deterministic(
        eccentric_pta_signal_planck18_1psr(
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
            log10_zc=log10_zc,
        ),
        name="gwecc",
    )

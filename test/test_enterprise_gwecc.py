import numpy as np
import pytest

from enterprise_gwecc import (
    eccentric_pta_signal_planck18,
    eccentric_pta_signal_planck18_1psr,
)

year = 365.25 * 24 * 3600
toas = np.linspace(0, 5 * year, 100)
pdist = 1.5
theta = 0.5
phi = 3.0
cos_gwtheta = 0.1
gwphi = 2.9
alpha = 0.3
psi = 1.2
cos_inc = 0.5
log10_M = 8.0
eta = 0.2
log10_F = -8.0
e0 = 0.3
gamma0 = gammap = 0.0
l0 = lp = 0.0
tref = max(toas)
log10_zc = -2.0


@pytest.mark.parametrize("psrTerm", [True, False])
def test_eccentric_pta_signal_planck18_1psr(psrTerm):
    res = eccentric_pta_signal_planck18_1psr(
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
        psrTerm=psrTerm,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize("psrTerm", [True, False])
def test_eccentric_pta_signal_planck18(psrTerm):
    res = eccentric_pta_signal_planck18(
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
        psrTerm=psrTerm,
    )
    assert np.all(np.isfinite(res))

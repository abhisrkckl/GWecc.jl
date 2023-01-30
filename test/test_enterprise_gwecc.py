import numpy as np
import pytest
import pathlib

from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA

from enterprise_gwecc import (
    eccentric_pta_signal_planck18,
    eccentric_pta_signal_planck18_1psr,
    gwecc_block,
    gwecc_1psr_block,
)


@pytest.fixture()
def psr():
    testdatadir = pathlib.Path(__file__).resolve().parent / "testdata"
    par = str(testdatadir / "J1909-3744_NANOGrav_12yv4.gls.par")
    tim = str(testdatadir / "J1909-3744_NANOGrav_12yv4.tim")
    return Pulsar(par, tim)


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


@pytest.mark.parametrize(
    "psrTerm, spline", [(True, True), (True, False), (False, True), (False, False)]
)
def test_eccentric_pta_signal_planck18_1psr(psrTerm, spline):
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
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, spline", [(True, True), (True, False), (False, True), (False, False)]
)
def test_eccentric_pta_signal_planck18(psrTerm, spline):
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
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, spline", [(True, True), (True, False), (False, True), (False, False)]
)
def test_gwecc_block(psr, psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_block(tref=tref, psrTerm=psrTerm, spline=spline)
    model = tm + wn + wf

    pta = PTA([model(psr)])

    assert len(pta.param_names) == (13 if psrTerm else 11)

    x0 = [param.sample() for param in pta.params]
    assert np.all(np.isfinite(x0))
    assert np.isfinite(pta.get_lnlikelihood(x0))
    assert np.isfinite(pta.get_lnprior(x0))


@pytest.mark.parametrize(
    "psrTerm, spline", [(True, True), (True, False), (False, True), (False, False)]
)
def test_gwecc_1psr_block(psr, psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_1psr_block(tref=tref, psrTerm=psrTerm, spline=spline)
    model = tm + wn + wf

    pta = PTA([model(psr)])

    assert len(pta.param_names) == (12 if psrTerm else 10)

    x0 = [param.sample() for param in pta.params]
    assert np.all(np.isfinite(x0))
    assert np.isfinite(pta.get_lnlikelihood(x0))
    assert np.isfinite(pta.get_lnprior(x0))

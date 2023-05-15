import numpy as np
import pytest
import pathlib
import itertools as it

from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA

from enterprise_gwecc import (
    eccentric_pta_signal,
    eccentric_pta_signal_1psr,
    eccentric_pta_signal_target,
    gwecc_block,
    gwecc_1psr_block,
    gwecc_target_block,
    gwecc_prior,
    gwecc_target_prior,
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
gwdist = 100.0
alpha = 0.3
psi = 1.2
cos_inc = 0.5
log10_M = 8.0
eta = 0.2
sigma = 0.3
rho = 0.5
log10_F = -8.0
e0 = 0.3
gamma0 = gammap = 0.0
l0 = lp = 0.0
tref = max(toas)
log10_A = -9.0
deltap = 100.0


@pytest.mark.parametrize("psrTerm, spline", it.product([True, False], [True, False]))
def test_eccentric_pta_signal_1psr(psrTerm, spline):
    res = eccentric_pta_signal_1psr(
        toas=toas,
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
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize("psrTerm, spline", it.product([True, False], [True, False]))
def test_eccentric_pta_signal(psrTerm, tie_psrTerm, spline):
    res = eccentric_pta_signal(
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
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, spline",
    it.product([True, False], [True, False]),
)
def test_eccentric_pta_signal_target(psrTerm, tie_psrTerm, spline):
    res = eccentric_pta_signal_target(
        toas,
        theta,
        phi,
        pdist,
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
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, tie_psrTerm, spline",
    it.product([True, False], [True, False], [True, False]),
)
def test_gwecc_block(psr, psrTerm, tie_psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_block(tref=tref, psrTerm=psrTerm, tie_psrTerm=tie_psrTerm, spline=spline)
    model = tm + wn + wf

    pta = PTA([model(psr)])

    if psrTerm:
        assert len(pta.param_names) == (13 if tie_psrTerm else 11)
    else:
        assert len(pta.param_names) == 11

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))


@pytest.mark.parametrize("psrTerm, spline", it.product([True, False], [True, False]))
def test_gwecc_1psr_block(psr, psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_1psr_block(tref=tref, psrTerm=psrTerm, spline=spline)
    model = tm + wn + wf

    pta = PTA([model(psr)])

    assert len(pta.param_names) == (9 if psrTerm else 8)

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))


@pytest.mark.parametrize(
    "psrTerm, tie_psrTerm, spline",
    it.product([True, False], [True, False], [True, False]),
)
def test_gwecc_target_block(psr, psrTerm, tie_psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_target_block(
        tref=tref,
        cos_gwtheta=cos_gwtheta,
        gwphi=gwphi,
        gwdist=gwdist,
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    model = tm + wn + wf

    pta = PTA([model(psr)])

    if psrTerm:
        assert len(pta.param_names) == (10 if tie_psrTerm else 8)
    else:
        assert len(pta.param_names) == 8

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_target_prior(pta, gwdist, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))

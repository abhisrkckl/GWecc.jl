"""Tests for the enterprise_gwecc python interface."""

import numpy as np
import pytest
import pathlib
import itertools as it
import glob

from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA
from enterprise.signals.parameter import Uniform

from enterprise_gwecc import (
    eccentric_pta_signal,
    eccentric_pta_signal_1psr,
    eccentric_pta_signal_target,
    gwecc_block,
    gwecc_1psr_block,
    gwecc_target_block,
    gwecc_prior,
    gwecc_target_prior,
    PsrDistPrior,
)

testdatadir = pathlib.Path(__file__).resolve().parent / "testdata"
par = str(testdatadir / "J1909-3744_NANOGrav_12yv4.gls.par")
tim = str(testdatadir / "J1909-3744_NANOGrav_12yv4.tim")
psr = Pulsar(par, tim)


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
psrdist = 1.0


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
def test_eccentric_pta_signal(psrTerm, spline):
    res = eccentric_pta_signal(
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
        psrTerm=psrTerm,
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, spline",
    it.product([True, False], [True, False]),
)
def test_eccentric_pta_signal_target(psrTerm, spline):
    res = eccentric_pta_signal_target(
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
        psrTerm=psrTerm,
        spline=spline,
    )
    assert np.all(np.isfinite(res))


@pytest.mark.parametrize(
    "psrTerm, tie_psrTerm, spline",
    it.product([True, False], [True, False], [True, False]),
)
def test_gwecc_block(psrTerm, tie_psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_block(
        tref=tref,
        psrdist=Uniform(1, 2),
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    model = tm + wn + wf

    pta = PTA([model(psr)])

    # assert len(pta.param_names) == (14 if (psrTerm and not tie_psrTerm) else 11)
    assert len(pta.param_names) == 11 + psrTerm * (1 + 2 * (not tie_psrTerm))

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))

    # Test that the prior will be zero if validation fails.
    wf = gwecc_block(
        tref=tref,
        log10_F=Uniform(-9, -2)("gwecc_log10_F"),
        psrdist=Uniform(1, 2),
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    model = tm + wn + wf
    pta = PTA([model(psr)])
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    log10_F_idx = pta.param_names.index("gwecc_log10_F")
    x0[log10_F_idx] = -3
    assert lnprior_fn(x0) == -np.inf


@pytest.mark.parametrize("psrTerm, spline", it.product([True, False], [True, False]))
def test_gwecc_1psr_block(psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_1psr_block(tref=tref, psrTerm=psrTerm, spline=spline)
    model = tm + wn + wf

    pta = PTA([model(psr)])

    assert len(pta.param_names) == 8 + (1 * psrTerm)

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))

    # Test that the prior will be zero if validation fails.
    wf = gwecc_1psr_block(
        tref=tref,
        log10_F=Uniform(-9, -2)("gwecc_log10_F"),
        psrTerm=psrTerm,
        spline=spline,
    )
    model = tm + wn + wf
    pta = PTA([model(psr)])
    lnprior_fn = gwecc_prior(pta, tref, tref, name="gwecc")
    log10_F_idx = pta.param_names.index("gwecc_log10_F")
    x0[log10_F_idx] = -3
    assert lnprior_fn(x0) == -np.inf


@pytest.mark.parametrize(
    "psrTerm, tie_psrTerm, spline",
    it.product([True, False], [True, False], [True, False]),
)
def test_gwecc_target_block(psrTerm, tie_psrTerm, spline):
    tref = max(psr.toas)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf = gwecc_target_block(
        tref=tref,
        cos_gwtheta=cos_gwtheta,
        gwphi=gwphi,
        gwdist=gwdist,
        psrdist=Uniform(1, 2),
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    model = tm + wn + wf

    pta = PTA([model(psr)])

    assert len(pta.param_names) == 8 + psrTerm * (1 + 2 * (not tie_psrTerm))

    x0 = [param.sample() for param in pta.params]
    lnprior_fn = gwecc_target_prior(pta, gwdist, log10_F, tref, tref, name="gwecc")
    if np.isfinite(lnprior_fn(x0)):
        assert np.all(np.isfinite(x0))
        assert np.isfinite(pta.get_lnlikelihood(x0))
        assert np.isfinite(pta.get_lnprior(x0))
    
    # Test that the prior will be zero if validation fails.
    wf = gwecc_target_block(
        tref=tref,
        cos_gwtheta=cos_gwtheta,
        gwphi=gwphi,
        gwdist=gwdist,
        log10_F=Uniform(-9, -2)("gwecc_log10_F"),
        psrdist=Uniform(1, 2),
        psrTerm=psrTerm,
        tie_psrTerm=tie_psrTerm,
        spline=spline,
    )
    model = tm + wn + wf
    pta = PTA([model(psr)])
    lnprior_fn = gwecc_target_prior(pta, gwdist, log10_F, tref, tref, name="gwecc")
    log10_F_idx = pta.param_names.index("gwecc_log10_F")
    x0[log10_F_idx] = -3
    assert lnprior_fn(x0) == -np.inf

def test_psrdist_prior():
    parfiles = sorted(glob.glob(f"{testdatadir}/*.par"))
    timfiles = sorted(glob.glob(f"{testdatadir}/*.tim"))

    psrs = [Pulsar(p, t) for p, t in zip(parfiles, timfiles)]

    psrdist_info = {
        "J0340+4130": [1.7115, 0.34230000000000005, "DM"],
        "J0613-0200": [1.0570824524312896, 0.1273862574811491, "PX"],
        "J0636+5128": [0.7272727272727273, 0.12429752066115701, "PX"],
        "J1909-3744": [1.1695906432748537, 0.01504736500119695, "PX"],
    }

    tref = max(max(psr.toas) for psr in psrs)

    tm = MarginalizingTimingModel()
    wn = MeasurementNoise(efac=1.0)
    wf1 = gwecc_block(
        tref=tref,
        psrdist=PsrDistPrior(psrdist_info),
        psrTerm=True,
        tie_psrTerm=True,
        spline=False,
    )
    model = tm + wn + wf1

    pta = PTA([model(psr) for psr in psrs])
    assert len(pta.param_names) == 11 + len(psrs)

    x0 = [p.sample() for p in pta.params]
    assert all(np.isfinite(x0))
    assert np.isfinite(pta.get_lnprior(x0))

    wf2 = gwecc_target_block(
        tref=tref,
        cos_gwtheta=cos_gwtheta,
        gwphi=gwphi,
        gwdist=gwdist,
        psrdist=PsrDistPrior(psrdist_info),
        psrTerm=True,
        tie_psrTerm=True,
        spline=False,
    )
    model = tm + wn + wf2

    pta = PTA([model(psr) for psr in psrs])
    assert len(pta.param_names) == 8 + len(psrs)

    x0 = [p.sample() for p in pta.params]
    assert all(np.isfinite(x0))
    assert np.isfinite(pta.get_lnprior(x0))

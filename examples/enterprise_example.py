"""This example demonstrates how to create an ENTERPRISE Signal object from GWecc."""

import matplotlib.pyplot as plt
import numpy as np
from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.signal_base import PTA

from enterprise_gwecc import eccentric_pta_signal_1psr, gwecc_1psr_block

# Calling eccentric_pta_signal_planck18_1psr
year = 365.25 * 24 * 3600
toas = np.linspace(0, 5 * year, 100)
pdist = 1.5
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
log10_dl = -15
res = eccentric_pta_signal_1psr(
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
    log10_dl,
)
plt.plot(toas, res)
plt.show()


# Creating a Signal object with gwecc_1psr_block
psr = Pulsar(
    "data/JPSR00_simulate.par",
    "data/JPSR00_simulate.tim",
    timing_package="tempo2",
)
tm = MarginalizingTimingModel()
ecw = gwecc_1psr_block(tref=max(psr.toas))
model = tm + ecw
pta = PTA([model(psr)])
print(pta.param_names)

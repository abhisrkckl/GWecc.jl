from enterprise.pulsar import Pulsar
import numpy as np
from Fe_stat import FeStat

psr = Pulsar("data/JPSR00_simulate.par", "data/JPSR00_simulate.tim")
noise_dict = {"JPSR00_efac": 1, "JPSR00_log10_t2equad": -20}
fe = FeStat([psr], noise_dict)


gw_skyloc = (0.5, 1.4)
log10_M = 9
eta = 0.2
log10_F = -8
e0 = 0.2
l0 = 0.0
F = fe.compute_Fe(gw_skyloc, log10_M, eta, log10_F, e0, l0)
print(F)
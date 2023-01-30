import numpy as np
import json

from enterprise.pulsar import Pulsar
from enterprise.signals.parameter import Uniform
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA
from enterprise_gwecc import gwecc_1psr_block
from matplotlib import pyplot as plt
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

parfile = "data/JPSR00_simulate.par"
timfile = "data/JPSR00_simulate.tim"
psr = Pulsar(parfile, timfile)

true_params = json.load(open("data/true_gwecc_params.dat", "r"))

# Only log10_zc is allowed to here for demonstration.
name = "gwecc"
tref = max(psr.toas)
priors = {
    "alpha": true_params["alpha"],  # Uniform(0, 1)(f"{name}_alpha"),
    "psi": true_params["psi"],  # Uniform(0, np.pi)(f"{name}_psi"),
    "cos_inc": true_params["cos_inc"],  # Uniform(-1, 1)(f"{name}_cos_inc"),
    "log10_M": true_params["log10_M"],  # Uniform(6, 9)(f"{name}_log10_M"),
    "eta": true_params["eta"],  # Uniform(0, 0.25)(f"{name}_eta"),
    "log10_F": true_params["log10_F"],  # Uniform(-9, -7)(f"{name}_log10_F"),
    "e0": true_params["e0"],  # Uniform(0.01, 0.8)(f"{name}_e0"),
    "gamma0": true_params["gamma0"],  # Uniform(0, np.pi)(f"{name}_gamma0"),
    "gammap": 0.0,  # Uniform(0, np.pi),
    "l0": true_params["l0"],  # Uniform(0, 2 * np.pi)(f"{name}_l0"),
    "lp": 0.0,  # Uniform(0, 2 * np.pi),
    "tref": tref,
    "log10_zc": Uniform(-5, -3)(f"{name}_log10_zc"),
}

tm = MarginalizingTimingModel()
wn = MeasurementNoise(efac=1)
wf = gwecc_1psr_block(**priors)
model = tm + wn + wf

pta = PTA([model(psr)])
print(pta.param_names)

x0 = np.array([p.sample() for p in pta.params])
print("Log-likelihoood at", x0, "is", pta.get_lnlikelihood(x0))

ndim = len(x0)
cov = np.diag(np.ones(ndim) * 0.01**2)
Nsteps = 10000
x0 = np.hstack(x0)
sampler = ptmcmc(
    ndim, pta.get_lnlikelihood, pta.get_lnprior, cov, outDir="gwecc_sims/chains/"
)
sampler.sample(x0, Nsteps)

chain_file = "gwecc_sims/chains/chain_1.txt"
chain = np.loadtxt(chain_file)

plt.subplot(121)
param_name = pta.param_names[0][len("gwecc_") :]
plt.plot(chain[:, :-4][:, 0])
plt.axhline(true_params[param_name], c="k")
plt.ylabel(param_name)

plt.subplot(122)
burn = 50
burned_chain = chain[burn:, :-4][:, 0]
plt.hist(burned_chain, bins=16)
plt.axvline(true_params[param_name], c="k")
plt.xlabel(param_name)
plt.show()

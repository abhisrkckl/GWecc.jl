import numpy as np
import json
import corner

from enterprise.pulsar import Pulsar
from enterprise.signals.parameter import Uniform
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA
from enterprise_gwecc import gwecc_1psr_block
from matplotlib import pyplot as plt
from PTMCMCSampler.PTMCMCSampler import PTSampler as ptmcmc

parfile = "gwecc_sims/JPSR00_simulate_1.par"
timfile = "gwecc_sims/JPSR00_simulate_1.tim"
try:
    psr = Pulsar(parfile, timfile)
except FileNotFoundError:
    print("Simulated par and tim files not found. Run simulation_example.py to create them.")

true_params = json.load(open("gwecc_sims/true_gwecc_params.dat", "r"))

# Only log10_dl is allowed to here for demonstration.
name = "gwecc"
tref = max(psr.toas)
priors = {
    "sigma": Uniform(0, np.pi)(f"{name}_sigma"), # true_params["sigma"], 
    "rho": true_params["rho"],  # Uniform(-np.pi, np.pi)(f"{name}_rho"),
    "log10_M": true_params["log10_M"],  # Uniform(6, 9)(f"{name}_log10_M"),
    "eta": true_params["eta"],  # Uniform(0, 0.25)(f"{name}_eta"),
    "log10_F": Uniform(-9, -7)(f"{name}_log10_F"),  # true_params["log10_F"],
    "e0": Uniform(0.01, 0.8)(f"{name}_e0"),  # true_params["e0"],
    "l0": Uniform(-np.pi, np.pi)(f"{name}_l0"),  # true_params["l0"],
    "tref": tref,
    "log10_A": Uniform(-11, -5)(f"{name}_log10_A"),
    "deltap": true_params["deltap"],
}

tm = MarginalizingTimingModel()
wn = MeasurementNoise(efac=1)
wf = gwecc_1psr_block(**priors)
model = tm + wn + wf

pta = PTA([model(psr)])
print(pta.param_names)

x0 = np.array([p.sample() for p in pta.params])
print("Log-likelihood at", x0, "is", pta.get_lnlikelihood(x0))

ndim = len(x0)
cov = np.diag(np.ones(ndim) * 0.01**2)
Niter = 1100000
x0 = np.hstack(x0)
sampler = ptmcmc(
    ndim,
    pta.get_lnlikelihood,
    pta.get_lnprior,
    cov,
    outDir="gwecc_sims/chains/",
    resume=False,
)
# This sometimes fails if the acor package is installed, but works otherwise.
# I don't know why.
sampler.sample(x0, Niter)

chain_file = "gwecc_sims/chains/chain_1.txt"
chain = np.loadtxt(chain_file)
print(chain.shape)

burn = chain.shape[0] // 3
burned_chain = chain[burn:, :-4]

for i in range(ndim):
    plt.subplot(ndim, 1, i + 1)
    param_name = pta.param_names[i][len("gwecc_") :]
    plt.plot(burned_chain[:, i])
    plt.axhline(true_params[param_name], c="k")
    plt.ylabel(param_name)
plt.show()

truths = [true_params[par[len("gwecc_") :]] for par in pta.param_names]
corner.corner(burned_chain, labels=pta.param_names, truths=truths)
plt.show()

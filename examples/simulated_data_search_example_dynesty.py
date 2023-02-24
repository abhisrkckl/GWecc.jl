import numpy as np
import dynesty
import pickle
import json

from enterprise.pulsar import Pulsar
from enterprise.signals.parameter import Uniform
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise.signals.white_signals import MeasurementNoise
from enterprise.signals.signal_base import PTA
from enterprise_gwecc import gwecc_1psr_block
from dynesty import plotting as dyplot
from matplotlib import pyplot as plt

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
    "log10_dl": Uniform(13, 15)(f"{name}_log10_dl"),
}

tm = MarginalizingTimingModel()
wn = MeasurementNoise(efac=1)
wf = gwecc_1psr_block(**priors)
model = tm + wn + wf

pta = PTA([model(psr)])
print(pta.param_names)

x0 = np.array([p.sample() for p in pta.params])
print(x0, pta.get_lnlikelihood(x0))


def prior_transform_fn(pta):
    mins = np.array([param.prior.func_kwargs["pmin"] for param in pta.params])
    maxs = np.array([param.prior.func_kwargs["pmax"] for param in pta.params])
    spans = maxs - mins

    def prior_transform(cube):
        return spans * cube + mins

    return prior_transform


prior_transform = prior_transform_fn(pta)

ndim = len(x0)
sampler = dynesty.DynamicNestedSampler(pta.get_lnlikelihood, prior_transform, ndim)
sampler.run_nested()
results = sampler.results
with open("JPSR00_simulate_result.pkl", "wb") as pkl:
    pickle.dump(results, pkl)

dyplot.cornerplot(results)
plt.savefig("JPSR00_simulate_result.pdf")
plt.show()

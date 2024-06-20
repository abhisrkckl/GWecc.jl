"""Injecting the eccentric SMBHB signal into TOAs using libstempo."""

import json
import os

import libstempo as lst
import libstempo.plot as lstplot
import libstempo.toasim as toasim
import matplotlib.pyplot as plt
import numpy as np

import enterprise_gwecc as gwecc

data_dir = "data"
parfile = f"{data_dir}/JPSR00_simulate.par"
timfile = f"{data_dir}/JPSR00_simulate.tim"

output_dir = "gwecc_sims"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

psr = lst.tempopulsar(parfile=parfile, timfile=timfile)
print(psr.name)


def save_psr_sim(psr, savedir):
    print("Writing simulated data for", psr.name)
    psr.savepar(f"{savedir}/{psr.name}_simulate_1.par")
    psr.savetim(f"{savedir}/{psr.name}_simulate_1.tim")
    lst.purgetim(f"{savedir}/{psr.name}_simulate_1.tim")


day_to_s = 24 * 3600
tref = float(max(psr.toas()) * day_to_s)
gwecc_params = {
    "sigma": 0.3,
    "rho": 0.6,
    "log10_M": 8.9,
    "eta": 0.25,
    "log10_F": -8.0,
    "e0": 0.5,
    "l0": 0.0,
    "tref": tref,
    "log10_A": -7,
    "deltap": 100.0,
}
with open(f"{output_dir}/true_gwecc_params.dat", "w") as outfile:
    json.dump(gwecc_params, outfile, indent=4)


# def get_pdist(jname):
#     dfile = f"{enterprise.__path__[0]}/datafiles/pulsar_distances.json"
#     with open(dfile, "r") as fl:
#         pdict = json.load(fl)
#         return pdict[jname][0]


def add_gwecc_1psr(psr, gwecc_params, psrTerm=False):
    toas = (psr.toas() * day_to_s).astype(float)

    signal = (
        np.array(
            gwecc.eccentric_pta_signal_1psr(
                toas=toas,
                psrTerm=psrTerm,
                **gwecc_params,
            )
        )
        / day_to_s
    )

    psr.stoas[:] += signal

    return signal


toasim.make_ideal(psr)
toasim.add_efac(psr, 1)
signal = add_gwecc_1psr(psr, gwecc_params)
print("Simulated TOAs for", psr.name)

lstplot.plotres(psr, label="Residuals")
plt.plot(psr.toas(), signal * day_to_s * 1e6, c="k", label="Injected signal")
plt.title(psr.name)
plt.legend()
plt.show()

psr.fit()
save_psr_sim(psr, output_dir)

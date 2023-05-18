import enterprise
import json
import pathlib
import glob

from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import TimingModel
from enterprise.signals.signal_base import PTA
from enterprise_gwecc import gwecc_target_block
from pdist_prior import set_pdist_from_dict, DeltaPsrDistPrior

with open("pulsar_distances_12.5.json", "r") as pdistfile:
    pdist_info = json.load(pdistfile)

mdc_path = pathlib.Path(enterprise.__file__).parent / "datafiles/mdc_open1/"
parfiles = sorted(glob.glob(str(mdc_path / "*.par")))
timfiles = sorted(glob.glob(str(mdc_path / "*.tim")))

psrs = [Pulsar(p, t) for p, t in zip(parfiles, timfiles)]

set_pdist_from_dict(psrs, pdist_info)

tref = max(max(psr.toas) for psr in psrs)
wf = gwecc_target_block(
    tref,
    cos_gwtheta=2.1,
    gwphi=1.3,
    gwdist=10.0,
    psi=0,
    cos_inc=0.5,
    eta=0.25,
    log10_F=-8,
    e0=0.3,
    gamma0=0,
    l0=0,
    log10_A=-8,
    delta_pdist=DeltaPsrDistPrior(pdist_info),
    psrTerm=True,
    tie_psrTerm=True,
    spline=True
)
tm = TimingModel()
model = tm + wf

pta = PTA([model(psr) for psr in psrs])
print(pta.param_names)

x0 = [p.sample() for p in pta.params]
print(list(zip(pta.param_names, x0)))
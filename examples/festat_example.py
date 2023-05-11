from enterprise.pulsar import Pulsar
from enterprise.signals.gp_signals import MarginalizingTimingModel
from enterprise_extensions.blocks import red_noise_block
from enterprise.signals.signal_base import PTA
from enterprise.signals.white_signals import MeasurementNoise
import numpy as np
from festat import PTAInnerProduct

psr = Pulsar("data/JPSR00_simulate.par", "data/JPSR00_simulate.tim")

tm = MarginalizingTimingModel()
wn = MeasurementNoise(efac=True)
rn = red_noise_block()
model = tm + wn + rn

pta = PTA([model(psr)])

params = [param.sample() for param in pta.params]
ptadot = PTAInnerProduct(pta, params)

print("rKr = ", ptadot.rKr)

a = [np.ones_like(res) for res in pta.get_residuals()]
print("aKa = ", ptadot.inner_product_aKb(a, a))

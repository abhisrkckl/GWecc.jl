import numpy as np
from scipy.special import erf
from enterprise.signals.parameter import Parameter, TruncNormalPrior, TruncNormalSampler
import enterprise
import json


def PsrDistPrior(psrdist_info: dict, dmdist_broaden_factor=2):
    """Truncated normal distribution for pulsar distances (psrdist).
    This is a truncated Normal distribution with mean and variance
    taken from the distance measurements, and a lower cutoff at 0, so
    that the pulsar distance doesn't go negative.

    The pulsar name is obtained by parsing the parameter name.
    The pulsar distance info is obtained from the following sources
    in that order:
        1. The `psrdist_info` dictionary
        2. The `pulsar_distances.json` file available in the ENTERPRISE
           distribution
        3. The default value used in ENTERPRISE (1.0 +/- 0.2)
    """

    class PsrDistPrior(Parameter):
        _pulsar_distance_info = psrdist_info
        _dmdist_broaden_factor = dmdist_broaden_factor
        _typename = "PsrDistPrior"
        _size = None

        def __init__(self, name):
            super().__init__(name)

            self.psrname = name.split("_")[0]

            if self.psrname in self._pulsar_distance_info:
                (
                    self.pdist,
                    self.pdist_sigma,
                    self.pdist_type,
                ) = self._pulsar_distance_info[self.psrname]
                assert self.pdist_type in ["PX", "DM"]

                if self.pdist_type == "DM":
                    # Broader distribution to account for systematic uncertainties
                    # in DM distance measurement
                    self.pdist_sigma *= self._dmdist_broaden_factor
            else:
                # Try with ENTERPRISE
                pdist_file_default = (
                    f"{enterprise.__path__[0]}/datafiles/pulsar_distances.json"
                )
                with open(pdist_file_default, "r") as f:
                    pdist_info_default = json.load(f)

                if self.psrname in pdist_info_default:
                    self.pdist, self.pdist_sigma = pdist_info_default[self.psrname]
                else:
                    # ENTERPRISE default value
                    self.pdist, self.pdist_sigma = 1, 0.2

            # Precompute the norm because erf is slow
            self.pnorm = (
                np.sqrt(2 / np.pi)
                / self.pdist_sigma
                / (1 + erf(self.pdist / np.sqrt(2) / self.pdist_sigma))
            )

        def _prior(self, name):
            """Hack to make sure Parameter.__init__ doesn't throw an error"""
            pass

        def prior(self, x):
            return TruncNormalPrior(
                x, self.pdist, self.pdist_sigma, 0, np.inf, norm=self.pnorm
            )

        def sample(self):
            return TruncNormalSampler(
                self.pdist, self.pdist_sigma, 0, np.inf, norm=self.pnorm
            )

        def __repr__(self):
            return self._typename

        @property
        def params(self):
            """Hack to make sure ENTERPRISE doesn't throw an error"""
            return [self]

    return PsrDistPrior

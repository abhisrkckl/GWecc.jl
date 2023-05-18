import numpy as np
from scipy.special import erf
from enterprise.signals.parameter import Parameter, TruncNormalPrior, TruncNormalSampler
import enterprise
import json

def set_pdist_from_dict(psrs: list, psrdist_info: dict, dmdist_broaden_factor=2):
    """Set pulsar distance information from a dictionary"""
    for psr in psrs:
        if psr.name in psrdist_info:
            Dp = psrdist_info[psr.name][0]
            Dp_sigma = psrdist_info[psr.name][1]
            Dp_type = psrdist_info[psr.name][2]
            assert Dp_type in ["DM", "PX"]

            if Dp_type == "DM":
                Dp_sigma *= dmdist_broaden_factor

            psr._pdist = (Dp, Dp_sigma)


# Interface
# pd = PsrDistPrior(psrdist_info) is a subclass of Parameter.
# pd(name) is a Parameter object. It figures out pulsar name from the parameter name.

# psrdist_info has entries like
# {
#    "JNAME" : (value, uncertainty, type)
# }


def DeltaPsrDistPrior(psrdist_info: dict, dmdist_broaden_factor=2):
    """Truncated normal distribution for pulsar distance correction (delta_pdist).
    This is a truncated Normal distribution with zero mean and unit variance. The
    lower cutoff is decided such that the pulsar distance does not go negative.
    """

    class DeltaPsrDistPrior(Parameter):
        _pulsar_distance_info = psrdist_info
        _dmdist_broaden_factor = dmdist_broaden_factor
        _typename = "PsrDistPrior"
        _size = None

        def __init__(self, name):
            super().__init__(name)

            self.psrname = name.split("_")[0]

            if self.psrname in self._pulsar_distance_info:
                self.pdist, self.pdist_sigma, self.pdist_type = self._pulsar_distance_info[
                    self.psrname
                ]
                assert self.pdist_type in ["PX", "DM"]

                if self.pdist_type == "DM":
                    self.pdist_sigma *= self._dmdist_broaden_factor
            else:
                # Try with ENTERPRISE
                pdist_file_default = enterprise.__path__[0] + "/datafiles/pulsar_distances.json"
                with open(pdist_file_default, "r") as f:
                    pdist_info_default = json.load(f)
                
                if self.psrname in pdist_info_default:
                    self.pdist, self.pdist_sigma = pdist_info_default[self.psrname]
                else:
                    # ENTERPRISE default value
                    self.pdist, self.pdist_sigma = 1, 0.2

            # Lower limit such that the pulsar distance is always positive.
            self.pmin = -self.pdist / self.pdist_sigma

            # Precompute the norm because erf is slow
            self.pnorm = np.sqrt(2 / np.pi) / (1 - erf(self.pmin / np.sqrt(2)))

        def _prior(self, name):
            """Hack to make sure Parameter.__init__ doesn't throw an error"""
            pass

        def prior(self, x):
            return TruncNormalPrior(x, 0, 1, self.pmin, np.inf, norm=self.pnorm)
        
        def sample(self):
            return TruncNormalSampler(0, 1, self.pmin, np.inf, norm=self.pnorm)

        def __repr__(self):
            return self._typename
        
        @property
        def params(self):
            """Hack to make sure ENTERPRISE doesn't throw an error"""
            return [self]

    return DeltaPsrDistPrior

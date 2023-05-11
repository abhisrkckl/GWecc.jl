from enterprise.signals.signal_base import PTA
from sksparse.cholmod import cholesky
import numpy as np
import scipy.sparse as sps
import scipy.linalg as sl


class PTAInnerProduct:
    """PTA inner product

    r = Timing residuals
    N = White noise matrix
    T = GP basis matrix
    phi = GP signal weights
    """

    def __init__(self, pta: PTA, params, noise_dict={}, cholesky_sparse=True):
        self.pta = pta
        self.noise_dict = noise_dict
        self.params = params if isinstance(params, dict) else pta.map_params(params)

        self.cholesky_sparse = cholesky_sparse

        self.pta.set_default_params(noise_dict)

        # This assertion makes sure that "r" is just the timing residuals.
        assert np.all(
            np.array(self.pta.get_delay(params)) == 0
        ), "The PTA object should not have any deterministic signals."

        self.Ns = self.pta.get_ndiag(params)
        self.Ts = self.pta.get_basis(params)
        self.rs = self.pta.get_residuals()

        self.TNTs = self.pta.get_TNT(params)
        self.TNrs = self.pta.get_TNr(params)
        self.phiinvs = self.pta.get_phiinv(params, logdet=False, method="cliques")
        self.rNrs = [rNr_logdet[0] for rNr_logdet in self.pta.get_rNr_logdet(params)]

        self.sum_rNr = np.sum(self.rNrs)

        self.rKr = self.inner_product_rKr()

    def inner_product_rKr(self):
        rKr = self.sum_rNr

        if self.pta._commonsignals:
            TNT = self.pta._lnlikelihood._block_TNT(self.TNTs)
            TNr = self.pta._lnlikelihood._block_TNr(self.TNrs)
            phiinv = self.phiinvs

            expval = (
                cholesky(TNT + sps.csc_matrix(phiinv))(TNr)
                if self.cholesky_sparse
                else sl.cho_solve(sl.cho_factor(TNT + phiinv), TNr)
            )
            rKr += np.dot(TNr, expval)
        else:
            for TNr, TNT, phiinv in zip(self.TNrs, self.TNTs, self.phiinvs):
                if TNr is None:
                    continue

                Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)

                expval = sl.cho_solve(sl.cho_factor(Sigma), TNr)

                rKr += np.dot(TNr, expval)

        return rKr

    def inner_product_aNb_arr(self, a, b):
        return [
            signalcollection.get_ndiag(self.params).solve(
                bb, left_array=aa, logdet=True
            )[0]
            for signalcollection, aa, bb in zip(self.pta._signalcollections, a, b)
        ]

    def inner_product_TNa_arr(self, a):
        return self.inner_product_aNb_arr(self.Ts, a)

    def inner_product_aKb(self, a, b):
        aNbs = self.inner_product_aNb_arr(a, b)
        TNas = self.inner_product_TNa_arr(a)
        TNbs = self.inner_product_TNa_arr(b)

        aKb = sum(aNbs)

        if self.pta._commonsignals:
            TNT = self.pta._lnlikelihood._block_TNT(self.TNTs)

            TNa = self.pta._lnlikelihood._block_TNr(TNas)
            TNb = self.pta._lnlikelihood._block_TNr(TNbs)

            phiinv = self.phiinvs

            expval = (
                cholesky(TNT + sps.csc_matrix(phiinv))(TNb)
                if self.cholesky_sparse
                else sl.cho_solve(sl.cho_factor(TNT + phiinv), TNb)
            )
            aKb += np.dot(TNa, expval)
        else:
            for TNa, TNb, TNT, phiinv in zip(TNas, TNbs, self.TNTs, self.phiinvs):
                Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)

                expval = sl.cho_solve(sl.cho_factor(Sigma), TNb)

                aKb += np.dot(TNa, expval)

        return aKb

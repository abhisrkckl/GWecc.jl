""" @author Amit Jit Singh, Nidhi Pant, Abhimanyu Susobhanan
"""

import numpy as np
import scipy.linalg as sl
from enterprise.signals.signal_base import PTA
from enterprise.signals.gp_signals import TimingModel
from enterprise.signals.parameter import Constant
from enterprise.signals.white_signals import MeasurementNoise
from enterprise_gwecc import eccentric_pta_signal_components

class FeStat:
    """
    Class for the Fe-statistic.

    :param psrs: List of `enterprise` Pulsar instances.
    :param noise_dict: Dictionary of noise parameters.

    """

    def __init__(self, psrs, noise_dict):
        self.psrs = psrs
        self.noise_dict = noise_dict

        efac = Constant()
        equad = Constant()
        
        ef = MeasurementNoise(efac=efac, log10_t2equad=equad)
        tm = TimingModel(use_svd=True)

        model = ef + tm

        self.pta = PTA([model(psr) for psr in psrs])

        # set white noise parameters
        self.pta.set_default_params(noise_dict)

        self.compute_Nmats()

    def compute_Nmats(self):
        """Makes the N matrix used in the F statistic"""
        self.TNTs = self.pta.get_TNT(self.noise_dict)
        self.phiinvs = self.pta.get_phiinv(
            self.noise_dict, logdet=False, method="partition"
        )
        # Get noise parameters for pta toaerr**2
        self.Nvecs = self.pta.get_ndiag(self.noise_dict)
        # Get the basis matrix
        self.Ts = self.pta.get_basis(self.noise_dict)

        self.Nmats = [
            make_Nmat(phiinv, TNT, Nvec, T)
            for phiinv, TNT, Nvec, T in zip(
                self.phiinvs, self.TNTs, self.Nvecs, self.Ts
            )
        ]

        self.Sigmas = [
            TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
            for TNT, phiinv in zip(self.TNTs, self.phiinvs)
        ]

        self.cfs = [sl.cho_factor(Sigma) for Sigma in self.Sigmas]

    def compute_Fe(
        self,
        gw_skyloc,
        log10_M,
        eta,
        log10_F,
        e0,
        l0,
        brave=False,
        maximized_parameters=False,
    ):
        """
        Computes the Fe-statistic (see Ellis, Siemens, Creighton 2012).

        :param log10_F: log10 GW frequency
        :param gw_skyloc: 2x{number of sky locations} array containing [theta, phi] for each queried sky location,
                          where theta=pi/2-DEC, phi=RA,
                          for single sky location use gw_skyloc= np.array([[theta,],[phi,]])
        :param brave: Skip sanity checks in linalg for speedup if True.
        :param maximized_parameters: Calculate maximized extrinsic parameters if True.

        :returns:
            fstat: value of the Fe-statistic
        :if maximized_parameters=True also returns:
            inc_max: Maximized value of inclination
            psi_max: Maximized value of polarization angle
            phase0_max: Maximized value of initial phase
            h_max: Maximized value of amplitude

        """
        psrs = self.psrs
        Ts = self.Ts
        Nmats = self.Nmats
        cfs = self.cfs

        npsr = len(self.psrs)
        N = np.zeros((npsr, 6))
        M = np.zeros((npsr, 6, 6))

        # fstat = np.zeros(gw_skyloc.shape[1])
        # if maximized_parameters:
        #     inc_max = np.zeros(gw_skyloc.shape[1])
        #     psi_max = np.zeros(gw_skyloc.shape[1])
        #     gamma0_max = np.zeros(gw_skyloc.shape[1])
        #     ampl_max = np.zeros(gw_skyloc.shape[1])

        cos_gwtheta, gwphi = np.cos(gw_skyloc[0]), gw_skyloc[1]
        tref = max(max(psr.toas) for psr in psrs)

        for idx, (psr, Nmat, cf, T) in enumerate(zip(psrs, Nmats, cfs, Ts)):

            # ntoa = len(psr.toas)

            # A = np.zeros((6, ntoa))
            A = np.array(eccentric_pta_signal_components(
                psr.toas,
                psr.theta,
                psr.phi,
                psr.pdist[0],
                cos_gwtheta,
                gwphi,
                log10_M,
                eta,
                log10_F,
                e0,
                l0,
                l0,
                tref,
                psrTerm=False
            ))

            # for i in range(6):
            #     A[i, :] = As[i]

            # ip1 = innerProduct_rr(A[0, :], psr.residuals, Nmat, T, cf, brave=brave)
            # ip2 = innerProduct_rr(A[1, :], psr.residuals, Nmat, T, cf, brave=brave)
            # ip3 = innerProduct_rr(A[2, :], psr.residuals, Nmat, T, cf, brave=brave)
            # ip4 = innerProduct_rr(A[3, :], psr.residuals, Nmat, T, cf, brave=brave)
            # ip5 = innerProduct_rr(A[4, :], psr.residuals, Nmat, T, cf, brave=brave)
            # ip6 = innerProduct_rr(A[5, :], psr.residuals, Nmat, T, cf, brave=brave)

            # N[idx, :] = np.array([ip1, ip2, ip3, ip4, ip5, ip6])
            N[idx, :] = [
                innerProduct_rr(A[jj, :], psr.residuals, Nmat, T, cf, brave=brave) 
                for jj in range(6)
            ]

            # define M matrix M_ij=(A_i|A_j)
            for jj in range(6):
                for kk in range(6):
                    M[idx, jj, kk] = innerProduct_rr(
                        A[jj, :], A[kk, :], Nmat, T, cf, brave=brave
                    )

        NN = np.copy(N)
        MM = np.copy(M)

        N_sum = np.sum(NN, axis=0)
        M_sum = np.sum(MM, axis=0)

        # take inverse of M
        Minv = np.linalg.pinv(M_sum)

        fstat = 0.5 * np.dot(N_sum, np.dot(Minv, N_sum))

        if maximized_parameters:
            a_hat = np.dot(Minv, N_sum)

            A_p = np.sqrt(
                (a_hat[0] + a_hat[4]) ** 2 + (a_hat[1] - a_hat[3]) ** 2
            ) + np.sqrt((a_hat[0] - a_hat[4]) ** 2 + (a_hat[1] + a_hat[3]) ** 2)
            A_c = np.sqrt(
                (a_hat[0] + a_hat[4]) ** 2 + (a_hat[1] - a_hat[3]) ** 2
            ) - np.sqrt((a_hat[0] - a_hat[4]) ** 2 + (a_hat[1] + a_hat[3]) ** 2)
            AA = A_p - np.sqrt(A_p**2 - A_c**2)
            # AA = A_p + np.sqrt(A_p**2 + A_c**2)

            # inc_max[j] = np.arccos(-A_c/AA)
            inc_max = np.arccos(A_c / AA)

            two_psi_max = np.arctan2(
                (A_c * a_hat[0] - A_p * a_hat[4]), (A_c * a_hat[3] + A_p * a_hat[1])
            )

            psi_max = 0.5 * np.arctan2(np.sin(two_psi_max), -np.cos(two_psi_max))

            # convert from [-pi, pi] convention to [0,2*pi] convention
            if psi_max < 0:
                psi_max += np.pi

            # correcting weird problem of degeneracy (psi-->pi-psi/2 and phi0-->2pi-phi0 keep everything the same)
            if psi_max > np.pi / 2:
                psi_max += -np.pi / 2

            half_phase0 = -0.5 * np.arctan2(
                A_p * a_hat[0] - A_c * a_hat[4], A_p * a_hat[3] + A_c * a_hat[1]
            )

            phase0_max = np.arctan2(
                -np.sin(2 * half_phase0), np.cos(2 * half_phase0)
            )

            # convert from [-pi, pi] convention to [0,2*pi] convention
            if phase0_max < 0:
                phase0_max += 2 * np.pi

        if maximized_parameters:
            return fstat, inc_max, psi_max, phase0_max
        else:
            return fstat


def innerProduct_rr(x, y, Nmat, Tmat, cf, TNx=None, TNy=None, brave=False):
    r"""
    Compute inner product using rank-reduced
    approximations for red noise/jitter
    Compute: x^T N^{-1} y - x^T N^{-1} T \Sigma^{-1} T^T N^{-1} y

    :param x: vector timeseries 1
    :param y: vector timeseries 2
    :param Nmat: white noise matrix
    :param Tmat: Modified design matrix including red noise/jitter
    :param Sigma: Sigma matrix (\varphi^{-1} + T^T N^{-1} T)
    :param TNx: T^T N^{-1} x precomputed
    :param TNy: T^T N^{-1} y precomputed
    :return: inner product (x|y)
    """

    # white noise term
    Ni = Nmat
    xNy = np.dot(np.dot(x, Ni), y)
    Nx, Ny = np.dot(Ni, x), np.dot(Ni, y)

    if TNx is None and TNy is None:
        TNx = np.dot(Tmat.T, Nx)
        TNy = np.dot(Tmat.T, Ny)

    if brave:
        # cf = sl.cho_factor(Sigma, check_finite=False)
        SigmaTNy = sl.cho_solve(cf, TNy, check_finite=False)
    else:
        # cf = sl.cho_factor(Sigma)
        SigmaTNy = sl.cho_solve(cf, TNy)

    ret = xNy - np.dot(TNx, SigmaTNy)

    return ret


def make_Nmat(phiinv, TNT, Nvec, T):

    Sigma = TNT + (np.diag(phiinv) if phiinv.ndim == 1 else phiinv)
    cf = sl.cho_factor(Sigma)
    # Nshape = np.shape(T)[0] # Not currently used in code

    TtN = np.multiply((1 / Nvec)[:, None], T).T

    # Put pulsar's autoerrors in a diagonal matrix
    Ndiag = np.diag(1 / Nvec)

    expval2 = sl.cho_solve(cf, TtN)
    # TtNt = np.transpose(TtN) # Not currently used in code

    # An Ntoa by Ntoa noise matrix to be used in expand dense matrix calculations earlier
    return Ndiag - np.dot(TtN.T, expval2)

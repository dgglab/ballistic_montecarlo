import numpy as np
import scipy.constants


class Bandstructure:
    def __init__(self, k, phi, B):
        self.k = k
        self.phi = phi

        self.rotate()
        self._k_to_r(B)

    def rotate(self):
        '''
        Rotates coordinate pairs by phi (positive is CCW)
        '''
        kx = np.cos(self.phi)*self.k[0] - np.sin(self.phi)*self.k[1]
        ky = np.cos(self.phi)*self.k[1] + np.sin(self.phi)*self.k[0]
        self.k = [kx, ky]

    def _k_to_r(self, B):
        '''
        Converts the Fermi surface to a real space trajectory according to B
        Input field in units of Tesla
        '''
        B_scaled = B*(scipy.constants.e) * \
            1E-16/(scipy.constants.hbar)  # T*C/hbar

        kx, ky = self.k
        kx, ky = [kx[::np.sign(B)], ky[::np.sign(B)]]
        rx = (1/B_scaled) * (ky - ky[0])
        ry = -(1/B_scaled) * (kx + kx[0])

        self.k = [kx, ky]
        self.r = [rx - np.mean(rx), ry - np.mean(ry)]
        self.dr = np.diff(self.r)

    def calculate_injection_prob(self, normal):
        '''
        computes the injection probabilities for a given state for
        an edge according to its normal
        '''

        S = np.sqrt(self.dr[0]**2 + self.dr[1]**2)
        S_max = np.max(S)
        theta = np.arctan2(self.dr[1], self.dr[0])

        prob = np.cos(theta - normal) * S/S_max
        prob = [0 if p < 0 else p for p in prob]
        prob = prob/np.sum(prob)
        return prob, np.cumsum(prob)

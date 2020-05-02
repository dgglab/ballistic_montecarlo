import numpy as np
import scipy.constants


class Bandstructure:
    def __init__(self, k, phi, B):
        # TODO check that first and last k are the same, probably needs some rounding
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

        # Here we are computing the integral across ks to back out the real
        # space trajectory. In the end, you find the following equations

        kx, ky = self.k
        kx, ky = [kx[::np.sign(B)], ky[::np.sign(B)]]
        rx = (1/B_scaled) * (ky - ky[0])
        ry = -(1/B_scaled) * (kx + kx[0])

        self.k = [kx, ky]
        self.r = [rx - np.mean(rx), ry - np.mean(ry)]
        self.dr = np.diff(self.r)

    def get_accent_color(self):
        return 'gotem'

    def calculate_injection_probs(self, edgenorms):
        '''
        computes the injection probabilities for a given state for
        an edge according to its normal
        '''

        # TODO fix keying of dict to be indexing
        # TODO normalize distribution

        self.in_prob = {}
        self.cum_prob = {}

        S = np.sqrt(self.dr[0]**2 + self.dr[1]**2)
        S_max = np.max(S)
        theta = np.arctan2(self.dr[1], self.dr[0])

        for norm in set(edgenorms):
            prob = np.cos(theta - norm) * S/S_max
            prob = [0 if p < 0 else p for p in prob]
            prob = prob/np.sum(prob)
            self.in_prob[norm] = prob
            self.cum_prob[norm] = np.cumsum(self.in_prob[norm])

    def get_injection_index(self, edgenorm):

        # TODO handle multiple edgenorms? Corner collison
        # TODO add fractional index injection

        return np.argmax(self.cum_prob[edgenorm] > np.random.rand())
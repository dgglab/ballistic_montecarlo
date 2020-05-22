import numpy as np
import scipy.constants


class Bandstructure:
    def __init__(self, k, phi, field):
        '''
        Inputs:
            k: list of two arrays (kx and ky) that form the Fermi surface. The states on the
               Fermi surface should be separated by arc length corresponding to the local
               Fermi velocity
            phi: angle to rotate the Fermi surface by
            field: Magnetic field used to make k to r
        '''
        self.k = k
        self.phi = phi
        self.rotate()
        self._k_to_r(field)

    def rotate(self):
        '''
        Rotates coordinate pairs by phi (positive is CCW)
        '''
        kx = np.cos(self.phi)*self.k[0] - np.sin(self.phi)*self.k[1]
        ky = np.cos(self.phi)*self.k[1] + np.sin(self.phi)*self.k[0]
        self.k = [kx, ky]

    def _k_to_r(self, field):
        '''
        Converts the Fermi surface to a real space trajectory according to B
        Input field in units of Tesla
        '''
        field_scaled = field*(scipy.constants.e) * \
            1E-16/(scipy.constants.hbar)  # T*C/hbar

        kx, ky = self.k
        kx, ky = [kx[::int(np.sign(field))], ky[::int(np.sign(field))]]
        rx = (1/field_scaled) * (ky - ky[0])
        ry = -(1/field_scaled) * (kx + kx[0])

        self.k = [kx, ky]
        self.r = [rx - np.mean(rx), ry - np.mean(ry)]
        self.dr = np.diff(self.r)

    def calculate_injection_prob(self, normal_angle):
        '''
        computes the injection probabilities for a given state for
        an edge according to its normal angle
        '''

        S = np.sqrt(self.dr[0]**2 + self.dr[1]**2)
        S_max = np.max(S)
        theta = np.arctan2(self.dr[1], self.dr[0])

        prob = np.cos(theta - normal_angle) * S/S_max
        prob = [0 if p < 0 else p for p in prob]
        prob = prob/np.sum(prob)
        return prob, np.cumsum(prob)

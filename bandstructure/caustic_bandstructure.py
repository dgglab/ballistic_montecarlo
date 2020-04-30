import numpy as np


class Bandstructure:
    def __init__(self, k, phi):
        self.k = k
        self.phi = phi

        self.rotate()

    def rotate(self):
        # Rotates coordinate pairs by phi (positive is CCW)
        kx = np.cos(self.phi)*self.k[0] - np.sin(self.phi)*self.k[1]
        ky = np.cos(self.phi)*self.k[1] + np.sin(self.phi)*self.k[0]
        self.k = [kx, ky]

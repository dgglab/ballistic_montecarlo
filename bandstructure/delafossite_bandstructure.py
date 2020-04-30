import numpy as np


def delafossite(n=1000):
    '''
    Returns a discretized Fermi surface for PdCoO2.
    In PdCoO2, the Fermi velocity is approximately constant so the Fermi
    surface consists of points separated by constant arc length.

    Optional inputs:
    n: number of points to generate on the Fermi surface
    '''
    k0 = 0.95  # A^-1
    k6 = 0.05  # A^-1
    k12 = 0.006  # A^-1
    vf = 4.8  # eVA (7.3E5 m/s)

    thetas = np.linspace(0, 2*np.pi, n)
    kf = k0 + k6 * np.cos(6 * thetas) + k12 * np.cos(12 * thetas)
    k = [kf*np.cos(thetas), kf*np.sin(thetas)]

    arclen = np.sqrt(np.sum(np.diff(k)**2, axis=0))
    arclen = np.insert(arclen, 0, 0)
    arclen_interp = np.linspace(0, np.sum(arclen), n)
    kx = np.interp(arclen_interp, np.cumsum(arclen), k[0])
    ky = np.interp(arclen_interp, np.cumsum(arclen), k[1])

    return [kx, ky]

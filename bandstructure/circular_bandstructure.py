import numpy as np


def circular(n=1000, k0=0.95):
    '''
    Returns a discretized circular Fermi surface
    The Fermi velocity is assumed constant so the Fermi
    sface consists of points separated by constant arc length.

    Optional inputs:
    n: number of points to generate on the Fermi surface
    k0: radius of Fermi surface in A^-1
    '''
    thetas = np.linspace(0, 2*np.pi, n)
    k = [k0*np.cos(thetas), k0*np.sin(thetas)]

    arclen = np.sqrt(np.sum(np.diff(k)**2, axis=0))
    arclen = np.insert(arclen, 0, 0)
    arclen_interp = np.linspace(0, np.sum(arclen), n)
    kx = np.interp(arclen_interp, np.cumsum(arclen), k[0])
    ky = np.interp(arclen_interp, np.cumsum(arclen), k[1])

    return [kx, ky]

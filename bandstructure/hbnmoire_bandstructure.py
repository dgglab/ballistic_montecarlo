import numpy as np


def hbnmoire(filename):
    '''
    Load in a file and returns a discretized Fermi surface for hBN/Gr moire.
    '''
    with open(filename) as f:
        ks = np.loadtxt(f, delimiter=',')

    kx, ky = zip(*ks)
    return [np.array(kx), np.array(ky)]

import numpy as np

def matriz_rigidez_local(E, A, L):
    return E*A/L* np.mat([[1, 0, -1, 0],
                          [0, 0, 0, 0],
                          [-1, 0, 1, 0],
                          [0, 0, 0, 0],])

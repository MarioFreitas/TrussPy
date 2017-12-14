import numpy as np

def matriz_rotacao(θ):
    return np.mat([[ np.cos(θ), np.sin(θ), 0, 0],
                   [-np.sin(θ), np.cos(θ), 0, 0],
                   [0, 0,  np.cos(θ), np.sin(θ)],
                   [0, 0, -np.sin(θ), np.cos(θ)],])


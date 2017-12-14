from copy import deepcopy as copy
import numpy as np

def rigidez_01(K_global, dir_restringidas):
    K_global_01 = copy(K_global)
    for i in dir_restringidas:
        K_global_01[i-1, :] = 0
        K_global_01[:, i-1] = 0

    for i in dir_restringidas:
        K_global_01[i-1, i-1] = 1

    return K_global_01

def deslocamentos_locais(X_global, correlacoes, mat_rot):
    x_local = np.zeros((4, 1))
    for i in correlacoes:
        x_local[i[0] - 1] = X_global[i[1] -1]
    return mat_rot * x_local

def deslocamentos_locais_simplificado(X_global, correlacoes, mat_rot):
    x_local = np.zeros((4, 1))
    for i in range(4):
        x_local[i] = X_global[correlacoes[i] -1]
    return mat_rot * x_local
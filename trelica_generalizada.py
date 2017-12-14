"""
Nome do Script: Solução de Treliça Pelo MEF
Autor: Mario Raul Freitas
Descrição: Esse script é uma demostração de como resolver uma treliça 2D
           utilizando o MEF. 
           Será utilizada uma treliça de aço (E = 205 GPa),
           com seção transversal circular, de diâmetro D = 10 cm.
           A treliça possui duas barras, com as seguintes coordenadas (xi, yi, xf, yf)
           Barra 1: (0, 0, 4, 0)
           Barra 2: (0, 3, 4, 0)
           Os apoios estão localizados nos pontos (0, 0) e (0, 3).
"""

# Importação das bibliotecas
import numpy as np
from matplotlib import pyplot as plt
from MEF import *
np.set_printoptions(linewidth=100, precision=2)

# Dados do problema
E = 205e9
D = 0.1
A = np.pi*(D**2)/4
coords_barras = [(0, 0, 4, 0),
                (0, 3, 4, 0)]
correlacoes = [[(1, 1), (2, 2), (3, 3), (4, 4)],
               [(1, 5), (2, 6), (3, 3), (4, 4)]]
dir_livres_globais = [3, 4]
dir_rest_globais = [1, 2, 5, 6]
F = np.array([0, 0, 0, -10000, 0, 0])

# Obtenção das propriedades geométricas
L = [comprimento(i) for i in coords_barras]
θ = [angulo(i) for i in coords_barras]

# Montagem das matrizes de rigidez locais
k_locais = [matriz_rigidez_local(E, A, i) for i in L]
r_locais = [matriz_rotacao(i) for i in θ]
k_locais_rot = [j.T * i * j for i, j in zip(k_locais, r_locais)]

# Montagem da matriz de rigidez global
nos = len(dir_livres_globais) + len(dir_rest_globais)
K_global = np.mat(np.zeros((nos, nos)))
for i, j in zip(k_locais_rot, correlacoes):
    alocacao(i, K_global, j)

# Técnica de Solução (1s e 0s)
K_global_01 = rigidez_01(K_global, dir_rest_globais)

# Cálculo dos deslocamentos globais
X_global = np.linalg.solve(K_global_01, F)

# Obtençao dos deslocamentos locais
x_locais = [deslocamentos_locais(X_global, i, j) for i, j in zip(correlacoes, r_locais)]

# Cálculo dos esforços
esforcos = [i * j for i, j in zip(k_locais, x_locais)]

# Cálculo das Reações
reacoes = K_global * np.mat(X_global).reshape((nos, 1))

# Resultados
print(f'Deslocamentos Globais\n{X_global}')
for i in range(len(esforcos)):
    print(f'\nEsforços Elemento {i+1}\n{esforcos[i].reshape((1,4))}')
print(f'\nForças e Reações Globais\n{reacoes.reshape(1, nos)}')
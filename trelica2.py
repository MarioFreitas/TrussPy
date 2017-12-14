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
coords_barra1 = (0, 0, 0, 4)
coords_barra2 = (0, 0, 4, 4)
coords_barra3 = (0, 0, 4, 0)
coords_barra4 = (0, 4, 4, 4)
coords_barra5 = (4, 0, 4, 4)
correlacoes_1 = [(1, 1), (2, 2), (3, 5), (4, 6)]    # [(local, global), (local, gloabal)]
correlacoes_2 = [(1, 1), (2, 2), (3, 7), (4, 8)]
correlacoes_3 = [(1, 1), (2, 2), (3, 3), (4, 4)]
correlacoes_4 = [(1, 5), (2, 6), (3, 7), (4, 8)]
correlacoes_5 = [(1, 3), (2, 4), (3, 7), (4, 8)]
dir_livres_globais = [7, 8]
dir_rest_globais = [1, 2, 3, 4, 5, 6]
F = np.array([0, 0, 0, 0, 0, 0, -10e3, -20e3])

# Obtenção das propriedades geométricas
L1 = comprimento(coords_barra1)
L2 = comprimento(coords_barra2)
L3 = comprimento(coords_barra3)
L4 = comprimento(coords_barra4)
L5 = comprimento(coords_barra5)
θ1 = angulo(coords_barra1)
θ2 = angulo(coords_barra2)
θ3 = angulo(coords_barra3)
θ4 = angulo(coords_barra4)
θ5 = angulo(coords_barra5)

# Montagem das matrizes de rigidez locais
k1 = matriz_rigidez_local(E, A, L1)
r1 = matriz_rotacao(θ1)
k1_ = r1.T*k1*r1

k2 = matriz_rigidez_local(E, A, L2)
r2 = matriz_rotacao(θ2)
k2_ = r2.T*k2*r2

k3 = matriz_rigidez_local(E, A, L3)
r3 = matriz_rotacao(θ3)
k3_ = r3.T*k3*r3

k4 = matriz_rigidez_local(E, A, L4)
r4 = matriz_rotacao(θ4)
k4_ = r4.T*k4*r4

k5 = matriz_rigidez_local(E, A, L5)
r5 = matriz_rotacao(θ5)
k5_ = r5.T*k5*r5

# Montagem da matriz de rigidez global
K_global = np.mat(np.zeros((8, 8)))
alocacao(k1_, K_global, correlacoes_1)
alocacao(k2_, K_global, correlacoes_2)
alocacao(k3_, K_global, correlacoes_3)
alocacao(k4_, K_global, correlacoes_4)
alocacao(k5_, K_global, correlacoes_5)

# Técnica de Solução (1s e 0s)
K_global_01 = rigidez_01(K_global, dir_rest_globais)

# Cálculo dos deslocamentos globais
X_global = np.linalg.solve(K_global_01, F)

# Obtençao dos deslocamentos locais
x1 = deslocamentos_locais(X_global, correlacoes_1, r1)
x2 = deslocamentos_locais(X_global, correlacoes_2, r2)
x3 = deslocamentos_locais(X_global, correlacoes_3, r3)
x4 = deslocamentos_locais(X_global, correlacoes_4, r4)
x5 = deslocamentos_locais(X_global, correlacoes_5, r5)

# Cálculo dos esforços
f1 = k1 * x1
f2 = k2 * x2
f3 = k3 * x3
f4 = k4 * x4
f5 = k5 * x5

# Cálculo das Reações
reacoes = K_global * np.mat(X_global).reshape((8, 1))

# Resultados
print(f'Deslocamentos Globais\n{X_global}')
print(f'\nEsforços Elemento 1\n{f1}')
print(f'\nEsforços Elemento 2\n{f2}')
print(f'\nEsforços Elemento 3\n{f3}')
print(f'\nEsforços Elemento 4\n{f4}')
print(f'\nEsforços Elemento 5\n{f5}')
print(f'\nForças e Reações Globais\n{reacoes}')
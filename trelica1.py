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
from MEF import *
np.set_printoptions(linewidth=100, precision=2)

# Dados do problema
E = 205e9
D = 0.1
A = np.pi*(D**2)/4
coords_barra1 = (0, 0, 4, 0)
coords_barra2 = (0, 3, 4, 0)
correlacoes_1 = [(1, 1), (2, 2), (3, 3), (4, 4)]    # [(local, global), (local, gloabal)]
correlacoes_2 = [(1, 5), (2, 6), (3, 3), (4, 4)]
dir_livres_globais = [3, 4]
dir_rest_globais = [1, 2, 5, 6]
F = np.array([0, 0, 0, -10000, 0, 0])

# Obtenção das propriedades geométricas
L1 = comprimento(coords_barra1)
L2 = comprimento(coords_barra2)
θ1 = angulo(coords_barra1)
θ2 = angulo(coords_barra2)

# Montagem das matrizes de rigidez locais
k1 = matriz_rigidez_local(E, A, L1)
r1 = matriz_rotacao(θ1)
k1_ = r1.T*k1*r1

k2 = matriz_rigidez_local(E, A, L2)
r2 = matriz_rotacao(θ2)
k2_ = r2.T*k2*r2

# Montagem da matriz de rigidez global
K_global = np.mat(np.zeros((6, 6)))
alocacao(k1_, K_global, correlacoes_1)
alocacao(k2_, K_global, correlacoes_2)

# Técnica de Solução (1s e 0s)
K_global_01 = rigidez_01(K_global, dir_rest_globais)

# Cálculo dos deslocamentos globais
X_global = np.linalg.solve(K_global_01, F)

# Obtençao dos deslocamentos locais
x1 = deslocamentos_locais(X_global, correlacoes_1, r1)
x2 = deslocamentos_locais(X_global, correlacoes_2, r2)

# Cálculo dos esforços
f1 = k1 * x1
f2 = k2 * x2

# Cálculo das Reações
reacoes = K_global * np.mat(X_global).reshape((6, 1))

print(reacoes, sep='\n\n')

# Resultados
# print(f'Deslocamentos Globais\n{X_global}')
# print(f'\nEsforços Elemento 1\n{f1}')
# print(f'\nEsforços Elemento 2\n{f2}')
# print(f'\nForças e Reações Globais\n{reacoes}')
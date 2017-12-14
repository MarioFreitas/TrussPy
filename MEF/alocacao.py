def alocacao(k_local, K_global, correlacoes):
    for i in correlacoes:
        for j in correlacoes:
            K_global[i[1]-1, j[1]-1] += k_local[i[0]-1, j[0]-1]


def alocacao_simplificada(k_local, K_global, correlacoes):
    for i in range(4):
        for j in range(4):
            K_global[correlacoes[i]-1, correlacoes[j]-1] += k_local[i, j]
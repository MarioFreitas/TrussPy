import numpy as np

def comprimento(coordenadas):
    xi = coordenadas[0]
    yi = coordenadas[1]
    xf = coordenadas[2]
    yf = coordenadas[3]

    return np.sqrt((xf - xi) ** 2 + (yf - yi) ** 2)

def angulo(coordenadas):
    xi = coordenadas[0]
    yi = coordenadas[1]
    xf = coordenadas[2]
    yf = coordenadas[3]

    x = xf - xi
    y = yf - yi

    if x == 0:
        return np.pi/2
    else:
        return np.arctan(y/x)
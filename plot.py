from matplotlib import pyplot as plt
import numpy as np
import json
from mpl_toolkits.mplot3d import Axes3D
from math import *


# аналитическое решение
def solve(x, y, z, t, Lx, Ly, Lz):
    at = pi * sqrt(1 / (Lx * Lx) + 1 / (Ly * Ly) + 1 / (Lz * Lz));
    return sin(pi / Lx * x) * sin(2 * pi / Ly * y) * sin(3 * pi / Lz * z) * cos(at * t + 2*pi);


def main():
    filename = 'layer.json'

    with open(filename) as f:
        data = json.load(f)

    Lx = data['Lx']
    Ly = data['Ly']
    Lz = data['Lz']
    t = data['t']
    N = data['N']

    points = data['points']

    x, y, z, u, u_analytical = [], [], [], [], []

    for xi, yi, zi, ui in data['points']:
        x.append(xi)
        y.append(yi)
        z.append(zi)
        u.append(ui)
        u_analytical.append(solve(xi, yi, zi, t, Lx, Ly, Lz))

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    img = ax.scatter(x, y, z, c=u, alpha=0.5, cmap='seismic')
    plt.savefig(f'numerical_N{N}_L{Lx}x{Ly}x{Lz}_t{t}.png', bbox_inches='tight')
    plt.clf()

    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
    img = ax.scatter(x, y, z, c=u_analytical, alpha=0.5, cmap='seismic')
    plt.savefig(f'analytical_N{N}_L{Lx}x{Ly}x{Lz}_t{t}.png', bbox_inches='tight')


if __name__ == '__main__':
    main()

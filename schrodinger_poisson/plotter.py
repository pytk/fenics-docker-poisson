from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import interp2d

class Constant(object):
    def __init__(self):
        self.Q = 1.60217662e-19
        self.HBAR = 6.582119514 * 10**-16
        self.H = 6.626070040 * 10**-34
        self.PI = 3.14159265359
        self.EPS = 8.854187817 * 10**-12
        self.T = 300
        self.M = 9.10938356 * 10**-31
        self.KB = 1.38064852* 10**-23

def electrostaticPotential(device, potential, iterate):
    X = np.linspace(device.xin, device.xfi, device.nx+1)
    Y = np.linspace(device.yin, device.yfi, device.ny+1)
    X, Y = np.meshgrid(X, Y)
    # plot wavefunction
    fig = plt.figure()
    ax = fig.gca(projection="3d")
    ax.plot_surface(X, Y, potential, linewidth=0.2, antialiased=True, cmap=plt.cm.coolwarm)
    ax.view_init(10, -120)
    plt.savefig("img/electrostatic_potential" + "_" + str(iterate) + ".png")

def electronDistribution(device, electron_distribution, iterate):
    fig = plt.figure()
    X = np.linspace(device.xin, device.xfi, device.nx+1)
    Y = np.linspace(device.yin, device.yfi, device.ny+1)
    X, Y = np.meshgrid(X, Y)
    plt.subplot(1,1,1)
    plt.pcolor(X, Y, electron_distribution)
    plt.savefig("img/electron_distribution" + "_" + str(iterate) + ".png")
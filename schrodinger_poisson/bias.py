import math

class Constant(object):
    def __init__(self):
        self.Q = 1.60217662e-19
        self.HBAR = 5.682119514e-16
        self.H = 4.135667662e-15
        self.PI = 3.14159265359
        self.EPS = 8.854187817 * 10**-21
        self.T = 300
        self.KB = 8.6173305* 10**-5
        self.M = 9.10938356 * 10**-31
        self.H = 4.135667662 * 10**-15

def bias(device, type, bias):
    """
    this is a function which return a potential for boundary condition
    in Non liner Poisson equation

    Args
    --------------
    device: Class
    type: string -> "Schottky" or "Ohmic"
    bias: applied voltage
    """

    constant = Constant()
    if (type == "Schottky"):
        potential = 0.0
        potential1 = bias

        potential2 = 0

        # schottky junction = workfunction of metal + conductionband minimum(band gap)
        # TaN = 4.05 [eV] workfunction
        potential3 = 4.05

        # germaium bandgap
        # Ge 0.66 [eV]

        potential3 -= 0.66
        potential3 *= -1

        potential = potential1 + potential2 + potential3
        return potential

    if (type == "Ohmic"):
        potential = 0
        potential1 = bias

        potential2 = constant.KB*constant.T/constant.Q

        nc = (constant.PI * constant.M * device.material["electron_effective_mass"] * constant.T / constant.H**2)
        nc = pow(nc, 1.5)

        # L vallye has 4 equivalent energy
        nc = nc * 2 * 4

        nv = (constant.PI * constant.M * device.material["hevy_hole_effective_mass"] * constant.T / constant.H**2)
        nv = pow(nv, 1.5)
        nv = 2 * nv

        ni = nc * nv * math.exp(-device.material["band_gap"]/constant.KB/constant.T)

        # this doner density should be reconsidered base on felmi level
        # this is a templary value
        doner = 1.0 * 10**3
        temp = math.sqrt(doner**2 + 4*ni**2)
        temp += doner
        temp /= 2*ni
        temp = math.log(temp)
        potential2 = potential2*temp

        # TaN = 4.05 [eV] workfunction
        potential3 = 4.05 

        potential = potential1 + potential2 + potential3
        return potential
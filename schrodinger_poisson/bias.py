from math import sqrt, exp, log, e

class Constant(object):
    def __init__(self):
        self.Q = 1.60217662 * 10**-19
        self.HBAR = 6.582119514 * 10**-16
        self.H = 6.626070040 * 10**-34
        self.PI = 3.14159265359
        self.EPS = 8.854187817 * 10**-12
        self.T = 300
        self.M = 9.10938356 * 10**-31
        self.KB = 1.38064852 * 10**-23

def fermiEnergy(material):
    """
    return fermi energy from flatt band energy
    """

    # TaN work function is 5.43 eV
    # Germanium work function is 4.05 eV
    # electron affinity is 1.2 eV
    # conduction band minimum is 0.66 eV (band gap)
    constant = Constant()
    electron_affinity = material["electron_affinity"]
    conduction_band_min = material["conduction_band_min"]
    metal_work_function = material["metal_work_function"]
    semiconductor_wave_function = material["workf_function"]
    flat_band = (metal_work_function - semiconductor_wave_function)
    fermi_energy = flat_band - metal_work_function + electron_affinity + conduction_band_min
    print("fermi energy : " + str(fermi_energy))
    return fermi_energy

def effectiveDensity(fermi_energy, material):
    """
    return effective carrier density from given fermi_energy
    """
    constant = Constant()
    conduction_band_min = material["conduction_band_min"]
    valence_band_min = material["valence_band_min"]
    effective_conduction_band_density = material["effective_conduction_band_density"]
    effective_valence_band_density = material["effective_valence_band_density"]
    intrinsic_carrier_concentration = material["intrinsic_carrier_concentration"]

    n_type = effective_conduction_band_density * exp((fermi_energy - conduction_band_min) / (constant.KB * constant.T / constant.Q))

    print("n-type = " + str(n_type))


    h_type = effective_valence_band_density * exp((valence_band_min - fermi_energy) / (constant.KB * constant.T /constant.Q))

    print("h-type = " + str(h_type))

    return n_type, h_type

"""
def thresholdVoltage(device, gate_voltage):
    # dielectric constant of SiO2 is 3.8 eV
    constant = Constant()
    material = device.material
    work_function = material["workf_function"]
    capacitance_oxyde = (constant.EPS * 3.8) / device.oxyde
    charge = -(gate_voltage - work_function) * capacitance_oxyde
    threshold_voltage = -charge / capacitance_oxyde + (0.5 * work_function)

    if threshold_voltage > gate_voltage:
        chargeA = -capacitance_oxyde*(gate_voltage - threshold_voltage)
        chargeB = - constant.Q * device.acceptor * device.oxyde
        threshold_voltage = -(chargeA + chargeB)/capacitance_oxyde + (0.5 * work_function)
        return threshold_voltage
    else:
        return threshold_voltage
"""
    


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
        potential3 = device.material["metal_work_function"]

        # germaium bandgap
        # Ge 0.66 [eV]

        potential3 -= device.material["band_gap"]
        potential3 *= -1

        potential = potential1 + potential2 + potential3
        print("schottky contact applied potential : " + str(potential))
        return potential

    if (type == "Ohmic"):
        # doner density
        doner = 1.0 * 10**21
        potential = 0
        potential1 = bias

        band_gap = device.material["band_gap"]

        potential2 = constant.KB*constant.T/constant.Q

        fermi_energy = fermiEnergy(device.material)

        n, p = effectiveDensity(fermi_energy, device.material)

        ni = sqrt(n * p) * exp(-band_gap/(2*constant.KB*constant.T/constant.Q))

        potential2 = potential2*log((sqrt(pow(doner, 2) * 4 * pow(ni, 2)) + doner) / (2 * ni))

        # TaN work function is 5.43 eV
        potential3 = -device.material["metal_work_function"]

        potential = potential1 + potential2 + potential3
        print("ohmic contact applied potential : " + str(potential))
        return potential
import numpy as np

def electricField(potential, device):
    """
    return electric field Ef from electrostatic potential
    Note: In this way to calculate electric field, this doesn't depend
    on wave function by schrodinger equation

    - electric field - dictionaly

    Args:
        - electrostatic potential ( poisson equation ) - numpy array
        - device (instance of Device Class which include device structure like length, width, etc ..)
    """

    # dimention of rectangle mesh
    electric_field = {}

    for index in range(device.ny + 1):
        x_potential = potential[index, :]

        gradient_potential = np.gradient(x_potential)

        electric_field[index] = gradient_potential

    return electric_field
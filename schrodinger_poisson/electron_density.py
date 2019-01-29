from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import scipy
import bias
import matplotlib.animation as animation
import mpmath
import constant
# Warning: from fenics import * will import both `sym` and
# `q` from FEniCS. We therefore import FEniCS first and then
# overwrite these objects.
from dolfin import *        

def subbandContribution(fermienergy, eigenvalue, cons, device):
    """
    to calculate number of electron in i-th subband dominated

    return 2d numpy array of electron number
    """
    hbar=1.05457266 * 10**-34
    constant_value = hbar*hbar*cons.PI
    effective_mass = device.material["electron_effective_mass"] * cons.M
    dos = effective_mass / constant_value
    state = dos*cons.KB*cons.T*math.log(math.exp((fermienergy + eigenvalue)/(cons.KB*cons.T)) + 1)
    return state

def electronDensity(eigenvector, state):
	"""
	print the electron density function(x)
	"""
	electron_density = []
	for i in range(len(eigenvector)):
		kth_term = [(j * j * state[i]) for j in eigenvector[i]]
		electron_density.append(kth_term)

	n = np.sum(electron_density, axis=0)
	return n

def effMassToArray(eigenvalue, potential, device):
	"""
	return an 1d numpy array for the structure giving the effective mass for a paticular energy.
	This methods implements Nelson's emprical 2-band model of non-parabolicity

	Args: 
		- eigenvalue : 1d numpy array of calculated eigen value from schrodinger.py
		- potential : 1d numpy array of calculated electrostatic potential from poisson.py
	"""
	effective_mass = device.material["electron_effective_mass"] * constant.M
	array = effective_mass * (1.0 + device.material["non-parabolicity"]*(eigenvalue-potential))
	return array
	

def calculationEffMass(wavefunction, eigenvalues , potential):
	"""
	find subband effective masses including non-parabolicy
	but stilling using a fixed effective mass for each subband dispersion
	"""
	conductionband_effmass_state = np.vstack([effMassToArray(eigenvalue, potential) for eigenvalue in eigenvalues])
	array = 1.0/np.sum(wavefunction**2/conductionband_effmass_state,axis=1)
	return array

def fermiLevel(density, eigenvalue, effective_mass_state):
	"""
	find the fermi level for 1d slice of density and eigen value
	"""
	hbar=1.05457266 * 10**-34
	energy_step = 1.0 * 10**-9

	def fermiLevel_0K(density, eigenvalue, effective_mass_state):
		Et, Ef = 0.0, 0.0
		for i, (Ei, csb_meff) in enumerate(zip(eigenvalue, effective_mass_state)):
			Et += Ei
			Efnew = (density*hbar**2*constant.PI/csb_meff+Et)/(i+1)
			if Efnew > Ei:
				Ef = Efnew
			else:
				break
				# found Ef and so we break

		number_of_state = np.zeros(eigenvalue.size)
		for i, (Ei, csb_meff) in enumerate(zip(eigenvalue, effective_mass_state)):
			temp = (Ef - Ei)*csb_meff/(hbar**2*pi)
			temp *= (Ni>0.0)
			number_of_state[i]=temp
		return Ef, number_of_state

	def fermiDirac(energy, fermienergy):
		"""
		return integral of Fermi Dirac Equation for energy independent density of states
		"""
		return constant.KB*constant.T*math.log(math.exp((fermienergy - energy)/(constant.KB*constant.T))+1)

	def function(fermienergy, effective_mass_state, density):
		"""
		return denisty
		"""
		temp = density
		for Ei, csb_meff in zip(eigenvalue, effective_mass_state):
			temp -= csb_meff*fermiDirac(Ei, fermienergy)/(hbar**2*constant.PI)
		return temp

	fermi_energy_0k, number_of_state_0k = fermiLevel_0K(density, eigenvalue, effective_mass_state)

	# implement Newton-Raphson method
	fermi_energy = fermi_energy_0k
	energy_step = energy_step

	while True:
		y = function(fermi_energy, effective_mass_state, density)
		dy = (function(fermi_energy+energy_step, effective_mass_state, density) - function(fermi_energy - energy_step, effective_mass_state, density))/(2.0 * energy_step)

	

	



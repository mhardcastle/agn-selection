# Code written by Google NotebookLM

import numpy as np
import astropy.units as u
from astropy.cosmology import LambdaCDM, z_at_value
from cosmology_change import *

def rho_l(L, z, rho_l0, L_lstar, alpha_l, z_l0, k_l):
  """
  Calculates the low-luminosity population density (Eq. 8 & 9).

  Args:
    L: Radio luminosity (W Hz^-1 sr^-1).
    z: Redshift.
    rho_l0: Normalization term for the low-luminosity population.
    L_lstar: Break luminosity for the low-luminosity population.
    alpha_l: Power-law slope for the low-luminosity population.
    z_l0: Maximum redshift for evolution of the low-luminosity population.
    k_l: Evolution parameter for the low-luminosity population.

  Returns:
    The low-luminosity population density.
  """
  if z < z_l0:
    return rho_l0 * (L / L_lstar)**(-alpha_l) * np.exp(-L / L_lstar) * (1 + z)**k_l
  else:
    return rho_l0 * (L / L_lstar)**(-alpha_l) * np.exp(-L / L_lstar) * (1 + z_l0)**k_l

def f_h(z, z_h0, z_h1, z_h2, model):
  """
  Calculates the high-luminosity evolution function (Eq. 11-13).

  Args:
    z: Redshift.
    z_h0: Peak redshift for the high-luminosity population.
    z_h1: Width of the Gaussian rise for the high-luminosity population.
    z_h2: Width of the Gaussian decline for the high-luminosity population (model C only).
    model: RLF model ('A', 'B', or 'C').

  Returns:
    The high-luminosity evolution function.
  """
  if model == 'A' or (model in ['B', 'C'] and z < z_h0):
    return np.exp(-0.5 * ((z - z_h0) / z_h1)**2)
  elif model == 'B' and z >= z_h0:
    return 1.0
  elif model == 'C' and z >= z_h0:
    return np.exp(-0.5 * ((z - z_h0) / z_h2)**2)
  else:
    raise ValueError("Invalid model specified.")

def rho_h(L, z, rho_h0, L_hstar, alpha_h, z_h0, z_h1, z_h2, model):
  """
  Calculates the high-luminosity population density (Eq. 10).

  Args:
    L: Radio luminosity (W Hz^-1 sr^-1).
    z: Redshift.
    rho_h0: Normalization term for the high-luminosity population.
    L_hstar: Break luminosity for the high-luminosity population.
    alpha_h: Power-law slope for the high-luminosity population.
    z_h0: Peak redshift for the high-luminosity population.
    z_h1: Width of the Gaussian rise for the high-luminosity population.
    z_h2: Width of the Gaussian decline for the high-luminosity population (model C only).
    model: RLF model ('A', 'B', or 'C').

  Returns:
    The high-luminosity population density.
  """
  return rho_h0 * (L / L_hstar)**(-alpha_h) * np.exp(-L_hstar / L) * f_h(z, z_h0, z_h1, z_h2, model)

def rho(L, z, rho_l0, L_lstar, alpha_l, z_l0, k_l, rho_h0, L_hstar, alpha_h, z_h0, z_h1, z_h2, model):
  """
  Calculates the total RLF (Eq. 7).

  Args:
    L: Radio luminosity (W Hz^-1 for H_0 = 70 km/s/Mpc)
    z: Redshift.
    rho_l0: Normalization term for the low-luminosity population.
    L_lstar: Break luminosity for the low-luminosity population.
    alpha_l: Power-law slope for the low-luminosity population.
    z_l0: Maximum redshift for evolution of the low-luminosity population.
    k_l: Evolution parameter for the low-luminosity population.
    rho_h0: Normalization term for the high-luminosity population.
    L_hstar: Break luminosity for the high-luminosity population.
    alpha_h: Power-law slope for the high-luminosity population.
    z_h0: Peak redshift for the high-luminosity population.
    z_h1: Width of the Gaussian rise for the high-luminosity population.
    z_h2: Width of the Gaussian decline for the high-luminosity population (model C only).
    model: RLF model ('A', 'B', or 'C').

  Returns:
    The total RLF.
  """
  # L here is for H_0 70 km/s/Mpc and in units of W/Hz
  # convert to 50 Mpc...
  L_scaled=convert_l(L,z,cosmo2,cosmo1)
  # convert to per sr
  L_scaled/=(4*np.pi)

  return 1.5*convert_rho(rho_l(L_scaled, z, rho_l0, L_lstar, alpha_l, z_l0, k_l) + rho_h(L_scaled, z, rho_h0, L_hstar, alpha_h, z_h0, z_h1, z_h2, model),z,cosmo1,cosmo2)

# Define parameters for Model C (Omega_M = 1) from Table 1
params_model_C_Omega1 = {
    'rho_l0': 10**(-7.120),
    'alpha_l': 0.539,
    'L_lstar': 10**(26.10),
    'z_l0': 0.706,
    'k_l': 4.30,
    'rho_h0': 10**(-6.196),
    'alpha_h': 2.27,
    'L_hstar': 10**(26.95),
    'z_h0': 1.91,
    'z_h1': 0.559,
    'z_h2': 1.378,
    'model': 'C'
}

params_model_C_Omega0 = {
    'rho_l0': 10**(-7.523),
    'alpha_l': 0.586,
    'L_lstar': 10**(26.48),
    'z_l0': 0.710,
    'k_l': 3.48,
    'rho_h0': 10**(-6.757),
    'alpha_h': 2.42,
    'L_hstar': 10**(27.39),
    'z_h0': 2.03,
    'z_h1': 0.568,
    'z_h2': 0.956,
    'model': 'C'
}

H0_1 = 50 * u.km / u.s / u.Mpc
H0_2 = 70 * u.km / u.s / u.Mpc

cosmo1 = LambdaCDM(H0=H0_1, Om0=0.0, Ode0=0)
cosmo2 = LambdaCDM(H0=H0_2, Om0=0.3, Ode0=0)

if __name__=='__main__':

  # Calculate RLF for a given luminosity and redshift
  L = np.logspace(21,29,100) # Example luminosity
  import matplotlib.pyplot as plt

  for z in [0.01,0.5,1,2,3]:

      rlf = rho(L, z, **params_model_C_Omega0)
      #print(f"RLF at L = {L:.2e} W Hz^-1 sr^-1 and z = {z:.2f}: {rlf:.2e}")

      plt.plot(np.log10(L),np.log10(rlf),label=str(z))

  plt.xlim(21,29)
  plt.legend(loc=0)
  plt.show()

                              

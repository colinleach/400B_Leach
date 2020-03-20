# A collection of miscellaneous functions that don't fit into any of the other classes

# import modules
import numpy as np
from numpy.linalg import norm
import scipy.optimize as so

import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from galaxy.galaxy import Galaxy
from galaxy.galaxies import Galaxies
from galaxy.centerofmass import CenterOfMass

G_val = 4.498768e-6 # units of kpc^3/Gyr^2/Msun

# Function to compute the dynamical mass, given the observed size and velocity 
# dispersion of a galaxy
# See in-class lab 5

def wolf_mass(sigma, Re):
    """ 
    Wolf mass estimator from Wolf+ 2010

    Args:
        sigma : 
            1D line of sight velocity dispersion in km/s
        Re : 
            2D radius enclosing half the stellar mass in pc

    Returns: estimate of the dynamical mass within the half light radius in Msun
    """

    return 4 / G_val * sigma**2 * Re / 1000

# Stellar to Halo Mass Relation

# Following the work of Moster et al. 2013 (MNRAS, 428, 3121)
# https://ui.adsabs.harvard.edu/abs/2013MNRAS.428.3121M/abstract

# `Equation 2:` $ \frac{m}{M} = 2N \left[ \left( \frac{M}{M_1} \right)^{-\beta} + 
#               \left(\frac{M}{M_1} \right)^{\gamma} \right]$ 
# $m$ = stellar mass, $M$ = halo mass
# `Equation 11:` log $M_1(z) = M_{10} + M_{11} \frac{z}{z+1} $ 
# `Equation 12:` $N(z) = N_{10} + N_{11} \frac{z}{z+1} $
# `Equation 13:`  $\beta(z) = \beta_{10} + \beta_{11} \frac{z}{z+1} $
# `Equation 14:` $\gamma(z) = \gamma_{10} + \gamma_{11} \frac{z}{z+1} $

# See in-class lab 5

class AbundanceMatching:
    
    def __init__(self, M, z):
        " input: Halo mass (Msun) and Redshift"
        
        #initializing the parameters:
        self.M = M # Halo Mass in Msun
        self.z = z  # Redshift
        
        
    def logM1(self):
        """eq. 11 of Moster 2013
        input: 
            redshift
        output: 
            M1, characteristic mass in log(Msun)
        """

        M10      = 11.59
        M11      = 1.195 
        return M10 + M11*(self.z/(1+self.z))  
    
    
    def N(self):
        """eq. 12 of Moster 2013
        input: 
            redshift
        output: 
            Normalization for eq. 2
        """

        N10      = 0.0351
        N11      = -0.0247
    
        return N10 + N11*(self.z/(1+self.z))
    
    
    def Beta(self):
        """eq. 13 of Moster 2013
        input: 
            redshift
        output: 
            power of the low mass slope"""

        beta10      = 1.376
        beta11      = -0.826
    
        return beta10 + beta11*(self.z/(1+self.z))
    
    def Gamma(self):
        """eq. 14 of Moster 2013
        input: 
            redshift
        output: 
            power of the high mass slope """

        gamma10      = 0.608
        gamma11      = 0.329
    
        return gamma10 + gamma11*(self.z/(1+self.z))
    
    
    def SHMratio(self):
        """ 
        eq. 2 of Moster + 2013
        Inputs: 
            halo mass M in solar masses (NOT in logspce)
            redshift
        Outputs: 
            Stellar mass to halo mass ratio
        """

        M1 = 10**self.logM1() # Converting characteristic mass to Msun from Log(Msun)
        A = (self.M/M1)**(-self.Beta())  # Low mass end
        B = (self.M/M1)**(self.Gamma())   # High mass end
        Norm = 2*self.N() # Normalization
    
        SHMratio = Norm*(A+B)**(-1)
    
        return SHMratio
    
    
 # Q1: add a function to the class that takes the SHM ratio and returns 
# The stellar mass 
    def StellarMass(self):
        """
        Using eq 2 of Moster+ 2013

        Returns the stellar mass (M_sun)
        """
        
        return self.M * self.SHMratio()

#  Function that returns Sersic Profile for an Elliptical System
# See in-class lab 6

def sersic(R, Re, n, Mtot):
    """
    Input
        R:
            radius (kpc)
        Re:
            half mass radius (kpc)
        n:
            sersic index
        Mtot:
            total stellar mass

    Returns
        Surface Brightness profile in Lsun/kpc^2
    """

    # We are assuming M/L = 1, so the luminosity is:
    L = Mtot
    
    # the effective surface brightness is
    # Ie = L/7.2/pi/Re**2
    Ie = L / (7.2 * np.pi * Re**2)
        
    return Ie * np.exp(-7.67 * ((R/Re)**(1.0/n) - 1.0))


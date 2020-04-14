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

# from galaxy.galaxy import Galaxy
# from galaxy.galaxies import Galaxies
# from galaxy.centerofmass import CenterOfMass

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
    
    
    # Method that takes the SHM ratio and returns the stellar mass 

    def StellarMass(self):
        """
        Using eq 2 of Moster+ 2013

        Returns the stellar mass (M_sun)
        """
        
        return self.M * self.SHMratio()
        

def sersic(R, Re, n, Mtot):
    """
    Function that returns Sersic Profile for an Elliptical System
    (See in-class lab 6)

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

def HernquistM(r, a=60*u.kpc, M_halo=1.97e12*u.M_sun):
    """
    Args:
        r (Quantity, units of kpc): distance from center
        a (Quantity, units of kpc): scale radius
        M_halo (Quantity, units of M_sun): total DM mass 
        
    Returns:
        Total DM mass enclosed within r (M_sun)
    """
    
    return np.round(M_halo * r**2 / (a + r)**2, 2)

def jacobi_radius(r, M_host, M_sat):
    """
    The Jacobi Radius for a satellite on a circular orbit about an extended host, 
    where the host is assumed to be well modeled as an isothermal sphere halo:

    R_j = r * (M_sat / 2 M_host(<r))}^(1/3)

    For MW/LMC, the Isothermal Sphere approximation is not a bad one within 50 kpc.

    In other contexts, can be called the Roche radius, Roche limit or Hill radius.

    Args:
        r: 
            distance between stellite and host (kpc)
        M_host: 
            host mass enclosed within r (M_sun)
        M_sat: 
            satellite mass (M_sun)
    
    returns: 
        Jacobi radius (kpc)
    """

    return r * (M_sat / (2*M_host))**(1/3)

def jacobi_mass(Rj, r, Mhost):
    """
    Function that returns min mass of a satellite given its observed size + distance 
    from a massive host: Msat = (Rj/r)**3 * 2 * Mhost

    Args:
        Rj: 
            Jacobi radius (approx as observed size) (kpc)
        r: 
            distance between stellite and host (kpc)
        Mhost: 
            mass enclosed within r (M_sun)
    
    returns: 
        Minimum mass Msat of a satellite given its current size (M_sun)
    """
    
    return (Rj/r)**3 * 2 * M_host


# Some functions to calculate useful 3D rotation matrices

def rotation_matrix_to_vector(old_axis, to_axis=None):
    """
    Args: 
        old_axis (3-vector)
            Vector to be brought into alignment with `to_axis` by rotation about the origin
        to_axis (3-vector)
            Angular momentum vector will be aligned to this (default z_hat)

    Returns: 
        3x3 rotation matrix

    Based on Rodrigues' rotation formula
    Ref: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
    
    Note that orientation in the plane perpendicular to 'to_axis' is arbitrary 
    """

    if to_axis is None:
        to_axis = np.array([0, 0, 1])
    else:
        to_axis /= norm(to_axis) # we need a unit vector

    old_axis /= norm(old_axis)

    # cross product between old_axis and new axis
    k_vec = np.cross(old_axis, to_axis) # 3-vector
    s_sq = np.sum(k_vec**2) # scalar, sin theta

    # dot product between old_axis and new axis 
    c = np.dot(old_axis, to_axis) # scalar, cos theta

    # rotation matrix, 3x3
    kx, ky, kz = k_vec
    K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
    R = np.eye(3) + K + K@K * (1 - c) / s_sq

    return R

def z_rotation_matrix(pt1, pt2):
    """
    Rotates about z-axis to line up two given points along the x-axis

    Args:
        pt1, pt2 (2-component iterables)
            define points to be placed on the  x-axis

    Returns:
        3x3 rotation matrix
    """

    diff = pt2 - pt1
    theta = -np.arctan(diff[1] / diff[0])
    print(theta*180/np.pi)
    R = np.array([[np.cos(theta), -np.sin(theta), 0], 
                  [np.sin(theta), np.cos(theta), 0],
                  [0, 0, 1]])
    return -R

def is_iterable(x):
    # a surprising omission from standard Python?
    try:
        iterator = iter(x)
    except TypeError:
        return False
    else:
        return True       

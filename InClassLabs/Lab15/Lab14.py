

# # In Class Lab 14 Solns
# # Cosmological Tools



# import modules
import numpy as np
import astropy.units as u

# For Lab 12: Import the constant for the speed of light
from astropy.constants import c
# For Lab 14: Import the boltzmann constant
from astropy.constants import k_B 

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib


# Lab 12 : added
# integrating tools from SciPy  
from scipy.integrate import simps








class CosmologicalTools:
    # Define a class that provides functions to compute various cosmological quantities
    # for a given cosmology  
        
    def __init__(self, OmegaM0, OmegaR0, OmegaL0, h):
        # initialize the instance of the class - for any given Cosmology
        # Input:    Omega M matter density parameter at z=0
        #           Omega R radiation density parameter at z=0
        #           Omega L  dark energy density parameter at z=0
        #           h  normalization for the hubble parameter at z=0
        
        # initialize the cosmology at z=0
        self.OmegaM0 = OmegaM0    ### Matter Density Parameter
        self.OmegaR0 = OmegaR0    ### Radiation Density Parameter
        self.OmegaL0 = OmegaL0    ### Dark Energy Density Parameter
        self.OmegaK0 = 1 - (OmegaM0 + OmegaR0 + OmegaL0)    #### Curvature Density Parameter
    
        self.h = h   # Normalization of Hubble Parameter   
        self.Ho = h*100*u.km/u.s/u.Mpc #  Hubble Constant at z=0  100 h km/s/Mpc
    
        self.k_B_ev = k_B.to(u.electronvolt/u.K)

    
    # Question 1 A)
    def HubbleParameter(self, z):
        # Function that defines the Hubble Parameter as a function of redshift
        # Input:   Redshift z 
        # Returns: The Hubble parameter at the given redshift in units of km/s/Mpc        
        
        # FILL THIS IN 
        M = self.OmegaM0*(1+z)**3
        R = self.OmegaR0*(1+z)**4
        L = self.OmegaL0
        K = self.OmegaK0*(1+z)**2
        
        return  self.Ho*np.sqrt(M+R+L+K)
    
    
    
    # Question 2 A)
    def OmegaM_Z(self,z):
        # Function that defines the matter density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Matter Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.OmegaM0*(1+z)**3*self.Ho**2/self.HubbleParameter(z)**2
    
    def OmegaR_Z(self,z):
        # Function that defines the radiation density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Radiation Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.OmegaR0*(1+z)**4*self.Ho**2/self.HubbleParameter(z)**2
    
    
    def OmegaL_Z(self,z):
        # Function that defines the dark energy density parameter as a function of redshift
        # Input:  Redshift z . Can be an array
        # Output:  Dark Energy Density Parameter at the given redshift.
        
        # FILL THIS IN
        return self.OmegaL0*self.Ho**2/self.HubbleParameter(z)**2
    
    
    # Question 1 A) 
    
    def LookBackTime(self, ze):
    # Function that computes the look back time at a given redshift
    # i.e. the difference in time from when a photon was emitted to when it is observed at present day.
    # Input:    Redshift emitted (ze). This cannot be an array. 
    # Output:   Time in units of Gyr Ago  (relative to present day). This is not an arrat    
    
        # Observed redshift  
        zo = 0
    
        # define an array with redshifts, spaced in intervals of 0.001 from zo to ze
        zrange = np.arange(zo, ze, 1e-4)
        
        # y = (1/H(zrange)).to(GYR)  /  (1+zrange)
        # But need to correct units of 1/H to be Gyr rather than seconds  
        # use the astropy.units functionality .to(units)
        # FILL THIS IN 
        y = (1.0/self.HubbleParameter(zrange)).to(u.Gyr)/(1+zrange)  
        
        # Integrate y numerically over zrange and return in units of Gyr
        # FILL THIS IN 
        return simps(y,zrange)*u.Gyr
    
    
    # Question 2 A) 
    
    def ComovingDistance(self, ze):
    # Function that returns the Comoving Radial Distance to an object at a given redshift
    # Distance to a galaxy that is moving with the Hubble Flow (expanding universe) at a given redshift
    # Input:    Redshift observed (zo) 
    #           Redshift emitted (ze)
    # Output:   DC in Mpc

        zo = 0
    
        # define an array with redshifts, spaced  in intervals of 0.001
        # Note that if you want redshifts smaller than 0.001 you'll need to refine this
        zrange = np.arange(zo,ze, 1e-4)
    
        # 1/H(zrange)*speed of light
        # Speed of light is loaded in modules from astropy, but in units of m/s --> need in km/s
        # FILL THIS IN
        y = c.to(u.km/u.s)*(1.0/self.HubbleParameter(zrange))
    
        # Integrate y numerically over zrange and return in units of Mpc
        # FILL THIS IN 
        return simps(y,zrange)*u.Mpc
    
    
    # Question 2 D) 
    
    def ProperDistance(self, zo, ze):
    # Function that returns the Proper Distance 
    # of an comoving distance measured today, at a given redshift (the distance measured by a ruler)
    # Input:    Redshift observed (zo) 
    #           Redshift of object (ze)
    # Output:   Proper Distance in Mpc
    
        # Comoving Distance (to emitted photon) [ independent of time] x the scale factor 
        # at the time of observation.
        return self.ComovingDistance(ze)/(1+zo)

 
    # Question 3 A)
    
    def LuminosityDistance(self, ze): 
    # Function that computes DL, the luminosity distance of a galaxy at a given redshift
    # Input:    Redshift emitted (ze) 
    # Output:   DL  in Mpc
        # this is an observable so
        zo = 0
        # Return  DL = DC*(1+z)
        return self.ComovingDistance(ze)*(1+ze)
    

    # Question 4 A)
    
    def AngularDiameterDistance(self, ze): 
    # Function that computes DA, the angular diameter distance at a given redshift
    # This is the proper distance between us and the source, at the time the photons were emitted.
    # Physical size of angular separation of 1 degree
    # Input:   Redshift emitted (ze)
    # Output:   DA  in Mpc
    
        
        # this is an observable so
        zo = 0
        
        # # FILL THIS IN
        # DA = DC/(1+z_emitted) = DL/(1+z)**2
        return self.ComovingDistance(ze)/(1+ze)     
    
    
    # Question 4 B) 
    
    def Separation(self, ze, angle):
    # Function to compute the physical distance corresponding to an angular separation at a given redshift
    # Input:    Redshift emmitted ze ,  
    #           angle: Angle between galaxies in arcsec
    # Output:  Distance in kpc
    
        # convert angle from arcsec to radians
        #    FILL THIS IN
        angleRad = (angle*u.arcsec).to(u.rad)
    
         # FILL THIS IN
        #   DA*angleRad
        return (self.AngularDiameterDistance(ze)*angleRad/u.rad).to(u.kpc)
    

 
    # Q1
    def Temperature(self, z):
        # function that returns the temperature of the universe as a function of redshift
        #Input: Redshift of interest
        # Returns: Temperature in K
        To = 2.73
        
        # Fill this in 
        # To(1+z)
        return To*(1+z)
    
    
    
    
    # Q2 
    # 
    ## Fill this in 
    def HorizonDistance(self, zo):
        # Compute the proper distance to the horizon at a given redshift
        # DC(zo,ze)/(1+zo)
        # input: Redshift of the observer
        # Returns the proper distance (Mpc)
        
        # define an array with redshifts, spaced  in intervals of 0.001
        # Note that if you want redshifts smaller than 0.001 you'll need to refine this
        zrange = np.arange(zo, 5000, 1e-3)
    
        # 1/H(zrange)*speed of light
        # Speed of light is loaded in modules from astropy, but in units of m/s --> need in km/s
        # FILL THIS IN
        y = c.to(u.km/u.s)*(1.0/self.HubbleParameter(zrange))
    
        # Integrate y numerically over zrange and return in units of Mpc
        # FILL THIS IN 
        Comoving = simps(y,zrange)*u.Mpc
        
        # Proper distance Comoving/(1+z)
        return Comoving/(1+zo)

    
    # Q3 
    ## Fill this in 
    def SoundHorizon(self,zdecouple):
        # The Maximal distance that sound can travel since the beginning of the universe to the time of 
        # decoupling
        # DC(zdecouple, zinfty)*speed of sound
        # Input:  Redshift of decoupling
        # Returns: Distance of Sound Horizon in Mpc
   
        
        # Proper distance HorizonDistance/sqrt(3)
        return self.HorizonDistance(zdecouple)/np.sqrt(3)
    







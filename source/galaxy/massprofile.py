import numpy as np
from numpy.linalg import norm

import astropy.units as u
from astropy.constants import G

from galaxy.centerofmass import CenterOfMass


class MassProfile:
    """
    Class to define mass enclosed as a function of radius and circular velocity
    profiles for a given galaxy and simulation snapshot

    Args:
        gal (Galaxy object):
            The desired galaxy/snap to operate on
        com_p (3-vector):
            Optional. The position of the disk CoM.
    """
    
    
    def __init__(self, gal, com_p=None):
        self.gal = gal
        
        self.com = CenterOfMass(gal)

        if com_p is None:
            self.com_p = self.com.com_p()
        else:
            self.com_p = com_p

        try:
            _ = self.com_p.unit
        except AttributeError: # not a Quantity
            self.com_p *= u.kpc
        
    def mass_enclosed(self, radii, ptype=None):
        """
        Calculate the mass within a given radius of the CoM 
        for a given type of particle.
        
        Args:
            radii (array of distances): spheres to integrate over
            ptype (int): particle type from (1,2,3), or None for total
            
        Returns:
            array of masses, in units of M_sun
        """
        
        if ptype is None:
            dataset = self.gal.data
        else:
            dataset = self.gal.filter_by_type(ptype)
            
        xyz = np.array([dataset[col] for col in ('x', 'y', 'z')]) * u.kpc
        
        # get coordinates centered on CoM
        xyz_centered = xyz - self.com_p[:, np.newaxis]
        
        # distances from CoM:
        dist = norm(xyz_centered, axis=0)
        
        within_r = lambda r: np.sum(dataset['m'][dist < r])
        masses = np.array([within_r(ri) for ri in radii]) * 1e10 * u.M_sun
        
        return masses
    
    def mass_enclosed_total(self, radii):
        """
        Calculate the mass within a given radius of the CoM, 
        summed for all types of particle.
        
        Args:
            radii (array of distances): spheres to integrate over
            
        Returns:
            array of masses, in units of M_sun
        """
        
#         return np.sum([self.mass_enclosed(ptype, radii) \
#                        for ptype in (1,2,3)], axis=0)
    
        return self.mass_enclosed(radii)
    
    def halo_mass(self):
        "Utility function to get a parameter for Hernquist mass"
        
        dataset = self.gal.filter_by_type(1) # just halo particles
        return np.sum(dataset['m']) * 1e10 * u.M_sun
    
    def hernquist_mass(self, r, a, M_halo=None):
        """
        Calculate the mass enclosed for a theoretical profile
        
        Args:
            r (Quantity, units of kpc): distance from center
            a (Quantity, units of kpc): scale radius
            M_halo (Quantity, units of M_sun): total DM mass (optional)

        Returns:
            Total DM mass enclosed within r (M_sun)
        """
        
        if M_halo is None:
            M_halo = self.halo_mass()

        return np.round(M_halo * r**2 / (a + r)**2, 2)    
    
    def circular_velocity(self, radii, ptype=None):
        """
        Calculate orbital velocity at a given radius from the CoM 
        for a given type of particle.
        
        Args:
            radii (array of distances): circular orbit
            ptype (int): particle type from (1,2,3), or None for total
            
        Returns:
            array of circular speeds, in units of km/s
        """
        
        if ptype is None:
            central_mass = self.mass_enclosed_total(radii)
        else:
            central_mass = self.mass_enclosed(radii, ptype)
            
        return np.sqrt(G * central_mass / radii).to(u.km / u.s)
    
    def circular_velocity_total(self, radii):
        "Syntactic sugar for circular_velocity(radii, ptype=None)"
        
        return self.circular_velocity(radii)
    
    def circular_velocity_hernquist(self, radii, a, M_halo=None):
        central_mass = self.hernquist_mass(radii, a, M_halo)
        
        return np.sqrt(G * central_mass / radii).to(u.km / u.s)
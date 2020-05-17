import numpy as np
from numpy.linalg import norm
from scipy.optimize import curve_fit

import astropy.units as u
from astropy.constants import G

from galaxy.centerofmass import CenterOfMass
from galaxy.utilities import find_nearest


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
            r (Quantity, units of kpc): 
                distance from center
            a (Quantity, units of kpc): 
                scale radius
            M_halo (Quantity, units of M_sun): 
                total DM mass (optional)

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
            radii (array of distances): 
                circular orbit
            ptype (int): 
                particle type from (1,2,3), or None for total
            
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
        """
        Theoretical V_circ assuming halo mass follows a Hernquist profile

        Args:
            radii (array of distances): 
                circular orbit
            a (Quantity, units of kpc): 
                scale radius
            M_halo (Quantity, units of M_sun): 
                total DM mass (optional)               
        """

        central_mass = self.hernquist_mass(radii, a, M_halo)
        
        return np.sqrt(G * central_mass / radii).to(u.km / u.s)

    def fit_hernquist_a(self, r_inner=1, r_outer=30, get_details=False):
        """
        Get `scipy.optimize` to do a non-linear least squares fit to find
        the best scale radius `a` for the Hernquist profile.

        Args:
            r_inner (numeric):
                Optional. Minimum radius to consider (implicit kpc). 
                Avoid values < 1 as they cause numeric problems.
            r_outer (numeric):
                Optional. Maximum radius to consider (implicit kpc). 
        """

        # A function suitable for curve_fit. Units must be removed.
        def hq(r, a):
            masses = (self.hernquist_mass(r*u.kpc, a*u.kpc)).value
            return np.log(masses)

        # The fitting has problems inside 1 kpc, so use a more restricted 
        # range of radii than for the plots
        radii_outer = np.linspace(r_inner, r_outer) * u.kpc

        # y values to fit to
        halo = np.log(self.mass_enclosed(radii_outer, 1).value)
        
        # run the fit and store the optimum a value
        popt, pcov = curve_fit(hq, radii_outer, halo, (60,))
        fitted_a = np.round(popt[0], 1)*u.kpc
        perr = np.round(np.sqrt(np.diag(pcov))[0], 1)*u.kpc

        if get_details:
            return fitted_a, perr
        else:
            return fitted_a

    def sersic(self, R, Re, n, Mtot):
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

    def bulge_Re(self, R):
        """
        Find the radius enclosing half the bulge mass.

        Args:
            R (array of Quantity):
                Radii to consider (kpc)

        Returns:
            Re (Quantity) :
                Radius enclosing half light/mass (kpc)
            bulge_total (numeric):
                Mass of entire bulge (M_sun, no units)
            bulgeI (array of Quantity):
                Surface brightness at radii R (kpc^-2), assuming M/L=1
        """

        bulge_mass = self.mass_enclosed(R, 3).value
        bulgeI = bulge_mass / (4 * np.pi * R**2) # divisor is area of sphere
        bulge_total = np.max(bulge_mass)
        Blow = bulge_total / 2.0
        # Bhigh = bulge_total / 2.0 + bulge_total / 2.0 * 0.2
        # index = np.where((bulge_mass > Blow) & (bulge_mass < Bhigh))
        index = np.where((bulge_mass > Blow))
        Re_bulge = R[index][0]
        
        return Re_bulge, bulge_total, bulgeI

    def fit_sersic_n(self, R, Re, bulge_total, bulgeI):
        """
        Get `scipy.optimize` to do a non-linear least squares fit to find
        the best value of `n` for a Sersic profile.

        Args:
            R (array of quantity):
                Radii at which to calculate fit (kpc)
            Re (Quantity) :
                Radius enclosing half light/mass (kpc)
            bulge_total (numeric):
                Mass of entire bulge (M_sun, no units)
            bulgeI (array of Quantity):
                Surface brightness at radii R (kpc^-2)

        Returns:
            best `n` value and error estimate
        """

        # A function suitable for curve_fit. Units must be removed.
        def ser(r, n):
            logI = self.sersic(r, Re.value, n, bulge_total)
            return np.log(logI)

        log_bulgeI = np.log(bulgeI.value)
        popt, pcov = curve_fit(ser, R, log_bulgeI, (60,))
        return popt[0], pcov[0][0]

    def density_profile_shell(self, radii, m, xyz):
        """
        Calculates mass density in successive spherical shells

        Arg:
            radii (array of float):
                boundary values beteen shells (implicit kpc, no units)
            m (shape (N,) array of float):
                particle masses (implicit Msun, no units)
            xyz ((3,N) array of float):
                particle cartesian coordinates

        Returns:
            r_annuli: geometric mean of boundaries 
                (array is one shorter than input radii)
            rho: densities (Msun/kpc^3)
                (same length as r_annuli)
        """

        r = norm(xyz, axis=0)
        enc_mask = r[:, np.newaxis] < np.asarray(radii).flatten()
        m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)
        m_annuli = np.diff(m_enc) * 1e10 * u.M_sun
        rho = 3/(4*np.pi) * m_annuli / (radii[1:]**3 - radii[:-1]**3)
        r_annuli = np.sqrt(radii[1:] * radii[:-1])
        return r_annuli, rho

    def density_profile_sphere(self, radii, m, xyz):
        """
        Calculates average mass density within successive spherical radii

        Arg:
            radii (array of float):
                boundary values beteen shells (implicit kpc, no units)
            m (shape (N,) array of float):
                particle masses (implicit Msun, no units)
            xyz ((3,N) array of float):
                particle cartesian coordinates

        Returns:
            rho: densities (Msun/kpc^3)
                (same length as radii)
        """

        r = norm(xyz, axis=0)
        enc_mask = r[:, np.newaxis] < np.asarray(radii).flatten()
        m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)
        rho = 3/(4*np.pi) * m_enc / (radii**3)
        return rho

    def virial_radius(self, r_min=20, r_max=1000, rho_c=None):
        """
        Calculates radius where DM density falls to 200x critical density
        for the universe.

        Args:
            r_min, r_max (floats)
                optional, limits for search (implicit kpc, no units)
            rho_c (float or Quantity)
                optional, critical density for chosen cosmology

        Returns:
            r_200 (float): virial radius, implicit kpc
        """

        if rho_c is None:
            rho_c = 127.35344 # Msun/kpc^3
        else:
            try:
                rho_c = rho_c.to(u.Msun/u.kpc**3).value
            except:
                pass # already has no units

        radii = np.linspace(r_min, r_max, 200)

        m = self.gal.data['m'] * 1e10 # in Msun
        DM = np.where(self.gal.data['type']==1)
        m_dm = m[DM]

        com = CenterOfMass(self.gal, ptype=None)
        xyz, _ = com.center_com()
        xyz_dm = (xyz.T[DM]).T
        rho_av = self.density_profile_sphere(radii, m_dm, xyz_dm)

        idx, nearest = find_nearest(rho_av, 200*rho_c)
        r_200 = radii[idx]
        return r_200

    def virial_mass(self, r_200=None, ptype=None):
        """
        Mass enclosed by the virial radius
        """

        if r_200 is None:
            r_200 = self.virial_radius()

        return self.mass_enclosed([r_200*u.kpc,], ptype=ptype)

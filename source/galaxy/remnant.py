# standard Python imports
from pathlib import Path

# scientific package imports
import numpy as np
from numpy.linalg import norm, eigh
from scipy.optimize import curve_fit
import astropy.units as u

from galaxy.galaxy import Galaxy
from galaxy.utilities import is_iterable


class Remnant(Galaxy):
    """
    A class to work with the MW-M31 post-merger remant.

    Args:
        snap (int):
            Snap number, equivalent to time elapsed. 
            Defaults to the last timepoint.
        datadir (str):
            Directory to search first for the required file. Optional, and a
            default list of locations will be searched.
        usesql (bool):
            If True, data will be taken from a PostgreSQL database instead of
            text files.
        stride (int):
            Optional. For stride=n, get every nth row in the table.
            Only valid with usesql=True.
        ptype (int or iterable of int):
            Particle type: 1, 2, 3 or a combination of these

    Class attributes:
        data (np.ndarray):
            type, mass, position_xyz, velocity_xyz for each particle
    """

    def __init__(self, snap=801, datadir=None, usesql=False, stride=1, ptype=(2,3)):
        "Initial setup. Currently it calls read_file(), but this may change."

        self.snap = snap
        self.ptype = ptype

        if usesql:
            self.read_db(stride)
        else:
            raise NotImplementedError

    def read_db(self, stride):
        """
        Get relevant data from a PostgreSQL database and format it to be 
        identical to that read from test files. 

        Ex-disk and ex-bulge particles are included, not DM particles.

        Args:
            stride (int):
                Optional. For stride=n, get every nth row in the table.

        Changes:
            `self.time`, `self.particle_count` and `self.data` are set.

        Returns: nothing
        """

        from galaxy.db import DB

        db = DB()
        cur = db.get_cursor()

        # set the elapsed time
        sql_t = f"SELECT time FROM simdata WHERE galname in ('MW', 'M31')"
        sql_t += f" and snap={self.snap} LIMIT 1"
        cur.execute(sql_t)
        time = cur.fetchone()
        try:
            self.time = time[0] * u.Myr
        except TypeError:
            print(self.name, self.snap, ptype)

        # set the bulk of the data
        colheads = ','.join(['galname','type','m','x','y','z','vx','vy','vz'])
        if stride > 1:
            sql_d = f"SELECT {colheads}, ROW_NUMBER() OVER () as rn from simdata"
        else:
            sql_d = f"SELECT {colheads} from simdata"
        sql_d += f"  where galname in ('MW', 'M31') and snap={self.snap}"
        if is_iterable(self.ptype):
            sql_d += f" and type in {self.ptype}"
        else:
            sql_d += f" and type={self.ptype}"
        sql_d += " ORDER BY galname, pnum"
        if stride > 1:
            sql_d = f"SELECT {colheads} from ( {sql_d} ) as t where rn % {stride} = 0" 

        dtype=[('galname', 'U3'), ('type', 'uint8'), ('m', '<f4'), 
                ('x', '<f4'), ('y', '<f4'), ('z', '<f4'), 
                ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        cur.execute(sql_d)
        self.data = np.array(cur.fetchall(), dtype=dtype)
        self.particle_count = len(self.data)
        
    def xyz(self):
        """
        Convenience method to get positions as a np.array of shape (3,N)
        """

        return np.array([self.data[xi] for xi in ('x','y','z')])

    def vxyz(self):
        """
        Convenience method to get velocities as a np.array of shape (3,N)
        """
        
        return np.array([self.data[vxi] for vxi in ('vx','vy','vz')])

    #
    # Methods to extimate ellipsoid principal axes
    #
    def I_tensor(self, m, x, y, z):
        """
        Args:
            m, x, y, z:
                1-D arrays with mass and coordinates (no units)
                
        Returns:
            3x3 array representing the moment of inertia tensor
        """
        
        # 3 moments of inertia for the diagonal
        Ixx = np.sum(m*(y**2 + z**2))
        Iyy = np.sum(m*(z**2 + x**2))
        Izz = np.sum(m*(x**2 + y**2))
        
        # 3 products of inertia for symmetric off-diagonals
        Ixy = Iyx = np.sum(m*x*y)
        Ixz = Izx = np.sum(m*x*z)
        Iyz = Izy = np.sum(m*y*z)

        # assemble the tensor and return it
        I = np.array([[Ixx, Ixy, Ixz],
                    [Iyx, Iyy, Iyz],
                    [Izx, Izy, Izz]])
        return I

    def ellipsoid_axes(self, m, x, y, z, r_lim=None):
        """
        Args:
            m, x, y, z:
                1-D arrays with mass and coordinates (no units)
            r_lim : float
                Radius to include in calculation (implicit kpc, no units)
                
        Returns:
            Two 3-tuples: relative semimajor axes and principal axis vectors
        """

        if r_lim is not None:
            r = np.sqrt(x**2 + y**2 + z**2)
            central = np.where(r < r_lim)
            x = x[central]
            y = y[central]
            z = z[central]
            m = m[central]

        # get moment-of-inertia tensor and eigenvalues/vectors
        I = self.I_tensor(m,x,y,z)
        w, v = eigh(I)

        # moments of intertia around principal axes:
        A, B, C = w / np.max(w)

        # principal axis units vectors are in columns of eigenvalues `v`
        vA, vB, vC = v.T

        # solve for semi-major axes of ellipsoid and normalize
        a = np.sqrt((B-A+C)/2)
        b = np.sqrt((C-B+A)/2)
        c = np.sqrt((A-C+B)/2)
        a, b, c = (a, b, c)/a

        return (a, b, c), (vA, vB, vC)

    #
    # Methods work with mass profiles, given CoM-centered coordinates
    #
    def sub_mass_enclosed(self, radii, m, xyz):
        """
        Calculate the mass within a given radius of the origin.
        Based on code in MassProfile, but this version assumes
        CoM-centric coordinates are supplied.

        Args:
            radii (array of distances): spheres to integrate over
            m (array of masses): shape (N,)
            xyz (array of Cartesian coordinates): shape (3,N)

        Returns:
            array of masses, in units of M_sun
        """

        # distances from CoM:
        dist = norm(xyz, axis=0) * u.kpc

        within_r = lambda r: np.sum(m[dist < r])
        masses = np.array([within_r(ri) for ri in radii]) * 1e10 * u.M_sun

        return masses   

    def hernquist_mass(self, r, a, M_halo):
        """
        Calculate the mass enclosed for a theoretical profile

        Args:
            r (Quantity, units of kpc): 
                distance from center
            a (Quantity, units of kpc): 
                scale radius
            M_halo (Quantity, units of M_sun): 
                total DM mass

        Returns:
            Total DM mass enclosed within r (M_sun)
        """

        return np.round(M_halo * r**2 / (a + r)**2, 2)    

    def fit_hernquist_a(self, m, xyz, r_inner=1, r_outer=100):
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
            masses = (self.hernquist_mass(r*u.kpc, a*u.kpc, M_halo)).value
            return np.log(masses)

        # The fitting has problems inside 1 kpc, so use a more restricted 
        # range of radii than for the plots
        radii_outer = np.linspace(r_inner, r_outer) * u.kpc
        
        M_halo = np.sum(m) * 1e10

        # y values to fit to
        halo = np.log(self.sub_mass_enclosed(radii_outer, m, xyz).value)

        # run the fit and store the optimum a value
        popt, pcov = curve_fit(hq, radii_outer, halo, (60,))
        fitted_a = np.round(popt[0], 3)*u.kpc
        perr = np.round(np.sqrt(np.diag(pcov))[0], 3)*u.kpc

        return fitted_a, perr

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

    def Re(self, R, m, xyz):
        """
        Find the radius enclosing half the mass.

        Args:
            R (array of Quantity):
                Radii to consider (kpc)
            m (array of float):
                masses (no units)
            xyz (array of float with shape (3,N)):
                coordinates (implicit kpc)
 
        Returns:
            sub_Re (Quantity) :
                Radius enclosing half light/mass (kpc)
            sub_total (numeric):
                Mass of entire system (M_sun, no units)
            subI (array of Quantity):
                Surface brightness at radii R (kpc^-2), assuming M/L=1
        """

        sub_mass = self.sub_mass_enclosed(R, m, xyz).value
        subI = sub_mass / (4 * np.pi * R**2)
        sub_total = np.max(sub_mass)
        Blow = sub_total / 2.0
        index = np.where((sub_mass > Blow))
        sub_Re = R[index][0]

        return sub_Re, sub_total, subI

    def fit_sersic_n(self, R, sub_Re, sub_total, subI):
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
            logI = self.sersic(r, sub_Re.value, n, sub_total)
            return np.log(logI)

        log_subI = np.log(subI.value)
        popt, pcov = curve_fit(ser, R, log_subI, (60,))
        return popt[0], pcov[0][0]


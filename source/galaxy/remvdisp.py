# import modules
import numpy as np
from numpy.linalg import norm
import astropy.units as u

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rcParams

from galaxy.remnant import Remnant
from galaxy.centerofmass import CenterOfMass
from galaxy.plots import Plots
from galaxy.timecourse import TimeCourse
from galaxy.utilities import rotation_matrix_to_vector


class Vdisp():
    """
    A class to work with rotation and velocity dispersion of the MW-M31 
    post-merger remant.

    Args:
        snap (int):
            Snap number, equivalent to time elapsed. 
            Defaults to the last timepoint.
        usesql (bool):
            If True, data will be taken from a PostgreSQL database instead of
            text files.
        ptype (int or iterable of int):
            Particle type: 1, 2, 3 or a combination of these
    """
    
    def __init__(self, snap=801, ptype=(2,3), r_lim=60):

        self.snap = snap
        tc = TimeCourse()
        self.t = tc.snap2time(snap)
        self.p = Plots()
        
        self.remnant = Remnant(snap=snap, ptype=ptype, usesql=True)
        self.calc_centered()
        self.rotate(r_lim)
        
        self.xbins = None
        self.means_yx = None
        self.means_zx = None

    def calc_centered(self):
        """
        Sets the CoM position and velocity, plus particle coordinates centered on the CoM.

        No args, no return value.
        """

        # use all stars, not just disk
        self.com = CenterOfMass(self.remnant, ptype=None) 
        self.com_p = self.com.com_p()
        self.com_v = self.com.com_v(self.com_p)
        self.pos, self.v = self.com.center_com()
        
    def rotate(self, r_lim):
        """
        Creates transformed coordinates with the angular momentum vector along z.

        Arg:
            r_lim (float):
                Only consider particles within this radius when computing L-hat. 
                Implicit kpc, no units.
        """

        L, _, _ = self.com.angular_momentum(self.com_p.value, self.com_v.value, r_lim=r_lim)
        R = rotation_matrix_to_vector(L)

        self.xyz_rot = R @ self.pos
        self.vxyz_rot = R @ self.v
        

    def calc_v_sigma(self, pos_index, v_index):
        """
        Calculate mean radial velocities and dispersions for bins along an axis.

        Args:
            pos_index, v_index (integers in (0,1,2)):
                Axis numbers for binning and for radial velocity.
                x=0, y=1, z=2

        Returns:
            Binned v_radial and dispersions, v_max and central dispersion. 
            All implicit km/s, no units.
        """

        means, sigmas = self.com.disp_by_radius(self.xyz_rot[pos_index], 
                                           self.vxyz_rot[v_index], self.xbins)
        vmax = (np.max(means) - np.min(means)) / 2
        sigma_central = np.max(sigmas)
        
#         print(f'snap: {snap}, vmax: {vmax:.2f}, sigma_central: {sigma_central:.2f},'
#               f' vmax/sigma: {vmax/sigma_central:.3f}')

        return means, sigmas, vmax, sigma_central

    def set_xbins(self, xbins):
        """
        Sets new x-boundaries for binning, invalidates any previous calculations.

        Args:
            xbins (array of float):
                Distances from CoM along chosen axis (signed, not just radius)
        """

        self.xbins = xbins
        self.means_yx = None
        self.means_zx = None
    
    def set_yx(self):
        "Convenience method to call calc_v_sigma() with y-axis and x-velocities"

        self.means_yx, self.sigmas_yx, self.vmax_yx, self.sigma_central_yx = self.calc_v_sigma(1, 0)
        
    def set_zx(self):
        "Convenience method to call calc_v_sigma() with z-axis and x-velocities"
        
        self.means_zx, self.sigmas_zx, self.vmax_zx, self.sigma_central_zx = self.calc_v_sigma(2, 0)
        
    
    def plot_yx(self, particles='Stellar', xlim=(-40,40), ylim1=(-120,120), ylim2=(0,200), 
                xbins=None, pngout=False, fname=None):
        "Wrapper for Plots.plot_v-sigma()"
    
        if self.means_yx is None:
            self.set_yx()
            
        title = f'{particles} ' + r'$\bar{v}$ and $\sigma$, y axis,' + f' t={self.t:.2f} Gyr'
        self.p.plot_v_sigma(self.xbins, self.means_yx, self.sigmas_yx, 
                       xlim=xlim, ylim1=ylim1, ylim2=ylim2, xlabel='y (kpc)', 
                       title=title, pngout=pngout, fname=fname)
        
    def plot_zx(self, particles='Stellar', xlim=(-40,40), ylim1=(-120,120), ylim2=(0,200), 
                xbins=None, pngout=False, fname=None):
        "Wrapper for Plots.plot_v-sigma()"

        if self.means_zx is None:
            self.set_zx()
            
        title = f'{particles} ' + r'$\bar{v}$ and $\sigma$, z axis,' + f' t={self.t:.2f} Gyr'
        self.p.plot_v_sigma(self.xbins, self.means_zx, self.sigmas_zx, 
                       xlim=xlim, ylim1=ylim1, ylim2=ylim2, xlabel='y (kpc)', 
                       title=title, pngout=pngout, fname=fname)
        
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
        """

        # use all stars, not just disk
        self.com = CenterOfMass(self.remnant, ptype=None) 
        self.com_p = self.com.com_p()
        self.com_v = self.com.com_v(self.com_p)
        self.pos, self.v = self.com.center_com()
        
    def rotate(self, r_lim):
        """
        """

        L, _, _ = self.com.angular_momentum(self.com_p.value, self.com_v.value, r_lim=r_lim)
        R = rotation_matrix_to_vector(L)

        self.xyz_rot = R @ self.pos
        self.vxyz_rot = R @ self.v
        

    def calc_v_sigma(self, pos_index, v_index):
        """
        """

        means, sigmas = self.com.disp_by_radius(self.xyz_rot[pos_index], 
                                           self.vxyz_rot[v_index], self.xbins)
        vmax = (np.max(means) - np.min(means)) / 2
        sigma_central = np.max(sigmas)
        
#         print(f'snap: {snap}, vmax: {vmax:.2f}, sigma_central: {sigma_central:.2f},'
#               f' vmax/sigma: {vmax/sigma_central:.3f}')

        return means, sigmas, vmax, sigma_central

    def set_xbins(self, xbins):
        self.xbins = xbins
        self.means_yx = None
        self.means_zx = None
    
    def set_yx(self):
        self.means_yx, self.sigmas_yx, self.vmax_yx, self.sigma_central_yx = self.calc_v_sigma(1, 0)
        
    def set_zx(self):
        self.means_zx, self.sigmas_zx, self.vmax_zx, self.sigma_central_zx = self.calc_v_sigma(2, 0)
        
    
    def plot_yx(self, particles='Stellar', xlim=(-40,40), ylim1=(-120,120), ylim2=(0,200), 
                xbins=None, pngout=False, fname=None):
    
        if self.means_yx is None:
            self.set_yx()
            
        title = f'{particles} ' + r'$\bar{v}$ and $\sigma$, y axis,' + f' t={self.t:.2f} Gyr'
        self.p.plot_v_sigma(self.xbins, self.means_yx, self.sigmas_yx, 
                       xlim=xlim, ylim1=ylim1, ylim2=ylim2, xlabel='y (kpc)', 
                       title=title, pngout=pngout, fname=fname)
        
    def plot_zx(self, particles='Stellar', xlim=(-40,40), ylim1=(-120,120), ylim2=(0,200), 
                xbins=None, pngout=False, fname=None):
    
        if self.means_zx is None:
            self.set_zx()
            
        title = f'{particles} ' + r'$\bar{v}$ and $\sigma$, z axis,' + f' t={self.t:.2f} Gyr'
        self.p.plot_v_sigma(self.xbins, self.means_zx, self.sigmas_zx, 
                       xlim=xlim, ylim1=ylim1, ylim2=ylim2, xlabel='y (kpc)', 
                       title=title, pngout=pngout, fname=fname)
        
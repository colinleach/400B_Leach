import numpy as np

import matplotlib.pyplot as plt
from matplotlib import rcParams

import astropy.units as u
import astropy.constants as c

from galaxy.galaxy import Galaxy
from galaxy.timecourse import TimeCourse
from galaxy.centerofmass import CenterOfMass


class SurfaceDensityProfile:
    """ 
    Calculate the surface density profile of a galaxy at a snapshot.
    Modified from code supplied by Rixin Li.
    """
    
    def __init__(self, gal, snap, radii=None, r_step=0.5, ptype=2, usesql=False):
        """ initialization
        input:
            galaxy (str): 
                the galaxy name: 'MW', 'M31', 'M33'
            snap_id (int): 
                the number of the snapshot of interest
            radii (iterable of float): 
                an array of radii from a near-center point to the outer skirt
                default value: np.arange(0.1, 0.95*r_outermost, r_step) kpc
            r_step (float):
                for defining a list of radii if none is supplied
        """
    
        galaxy = Galaxy(gal, snap, ptype=ptype, usesql=usesql)
        self.ptype = galaxy.data['type']
        
        tc = TimeCourse()
        self.t = tc.snap2time(snap) * u.Gyr
        self.gal = gal
        self.snap = snap

        self.com = CenterOfMass(galaxy, 2)
        self.com_r, self.com_v = tc.get_one_com(gal, snap) 

        # rotate the frame to align the disk with angular momentum
        self.alg_r, self.alg_v = self.com.rotate_frame(com_p=self.com_r, com_v=self.com_v)

        # calculate the radial distances and azimuthal angles in cylindrical coordinates
        self.cyl_r_mag = np.sqrt(np.sum(self.alg_r[:2,:]**2, axis=0))
        self.cyl_theta = np.arctan2(self.alg_r[1,:], self.alg_r[0,:]) * 180/np.pi

        # check if radii is already set
        if radii is None:
            self.radii = np.arange(0.1, 0.95 * self.cyl_r_mag.max(), r_step)
        else:
            self.radii = radii
            # can be improved by checking how many elements "radii" has

        # create the mask to select particles for each radius
        enc_mask = self.cyl_r_mag[:, np.newaxis] < np.asarray(self.radii).flatten()
        
        # calculate the enclosed masses within each radius
        # relevant particles will be selected by enc_mask (i.e., *1)
        # outer particles will be ignored (i.e., *0)
        self.m_enc = np.sum(self.com.m[:, np.newaxis] * enc_mask, axis=0)

        # use the difference between nearby elements to get mass in each annulus
        # N.B.: we ignored the very central tiny circle and a small portion of 
        # outermost particles
        self.m_annuli = np.diff(self.m_enc) # this array is one element less then m_enc
        
        # calculate the surface density by dividing the area of the annulus
        self.Sigma = self.m_annuli / (np.pi * (self.radii[1:]**2 - self.radii[:-1]**2))
        
        # we use the geometric mean of two consecutive elements in "radii" 
        # as the radius of each annulus
        # this array has the same number of elements as self.Sigma, 
        # can be used for plotting
        self.r_annuli = np.sqrt(self.radii[1:] * self.radii[:-1])

    def plot_xy(self, figsize=(8,8), xlim=(0,60), ylim=(0,60), 
                pfc='b', ngout=False, fname=None):
        """

        """

        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(self.alg_r[0,:], self.alg_r[1,:], ec="None", fc=fc, s=0.25)

        ax.set_xlabel('x (kpc)', fontsize=22)
        ax.set_ylabel("y (kpc)", fontsize=22)
        # ax.set_aspect=1.0
        ax.set_title(f"{self.gal} Stellar Disk (t={self.t:.2f})", fontsize=22)

        #set axis limits
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        #adjust tick label font size
        label_size = 22
        rcParams['xtick.labelsize'] = label_size 
        rcParams['ytick.labelsize'] = label_size

        fig.tight_layout()

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

    def plot_r_theta_scaled(self, figsize=(8,8), xlim=(0,60), ylim=(-180,180), 
                fc='b', pngout=False, fname=None):
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.scatter(self.cyl_r_mag, self.cyl_theta, ec="None", fc=fc, 
                s=0.1*np.log(self.cyl_r_mag**2))
        
        ax.set_xlabel('r (kpc)', fontsize=22)
        ax.set_ylabel(r"$\theta$ (deg)", fontsize=22)
        ax.set_title(f"{self.gal} Stellar Disk (t={self.t:.2f})", fontsize=22)

        #set axis limits
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        #adjust tick label font size
        label_size = 22
        rcParams['xtick.labelsize'] = label_size 
        rcParams['ytick.labelsize'] = label_size

        fig.tight_layout()

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

# import modules
import numpy as np
from numpy.linalg import norm
import scipy.optimize as so

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib import rcParams

class Plots():
    """
    """

    def __init__(self):
        pass

# Code for plotting contours
# from https://gist.github.com/adrn/3993992
# See in-class lab 7

    def density_contour(self, xdata, ydata, nbins_x, nbins_y, level_vals, 
                        ax=None, **contour_kwargs):
        """ 
        Create a density contour plot.
        
        Parameters
            xdata : numpy.ndarray

            ydata : numpy.ndarray

            nbins_x : int
                Number of bins along x dimension
            nbins_y : int
                Number of bins along y dimension
            level_vals : list of float
                Contour values to include
            ax : matplotlib.Axes (optional)
                If supplied, plot the contour to this axis. Otherwise, open a new figure
            contour_kwargs : dict
                kwargs to be passed to pyplot.contour()
                i.e. unknown number of keywords 
            
        Example Usage
            density_contour(x pos, y pos, contour res, contour res, level_vals, axis, 
            colors for contours)
            
            e.g.:

            colors=['red','orange', 'yellow']

            density_contour(xD, yD, 80, 80, [0.68, 0.8, 0.9], ax=ax, colors=colors)
        """

        def find_confidence_interval(x, pdf, confidence_level):
            return pdf[pdf > x].sum() - confidence_level

        H, xedges, yedges = np.histogram2d(xdata, ydata, bins=(nbins_x,nbins_y), normed=True)
        x_bin_sizes = (xedges[1:] - xedges[:-1]).reshape((1,nbins_x))
        y_bin_sizes = (yedges[1:] - yedges[:-1]).reshape((nbins_y,1))

        pdf = (H*(x_bin_sizes*y_bin_sizes))  
        
        X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
        Z = pdf.T  # transpose of distribution fxn
        fmt = {}
        
        ### Adjust Here #### 
        
        # Contour Levels Definitions 
        # brentq is root finding method
        lev_fn = lambda x : so.brentq(find_confidence_interval, 0., 1., args=(pdf, x))
        levels = [lev_fn(level) for level in level_vals[::-1]]
        
        # contour level labels
        strs = [str(level) for level in level_vals[::-1]]

        
        ###### 
        
        if ax == None:
            contour = plt.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
            for l, s in zip(contour.levels, strs):
                fmt[l] = s
            plt.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)

        else:
            contour = ax.contour(X, Y, Z, levels=levels, origin="lower", **contour_kwargs)
            for l, s in zip(contour.levels, strs):
                fmt[l] = s
            ax.clabel(contour, contour.levels, inline=True, fmt=fmt, fontsize=12)
        
        return contour

    def plot_density(self, rn, galname, snap, t, lim=None, nbins=150, 
                    pngout=False, fname=None):
        """
        Particle column density plot, 2D histograms overlaid with contours. 
        
        Shows x,y (face-on) and x,z (edge-on) orientations, assuming rotation 
        axis is along z.

        Args:
            rn (np.array, shape (3,N))
                (x,y,z) coordinates
            galname (str)
                Galaxy name, just for the title
            snap (int):
                Timepoint sequence number, as used in the data filename
            t (float):
                Time in Gyr
            lim (number):
                Optional, x and y limits of the plot in kpc
            nbins (int)
                Optional, for passing to plt.hist2d
            pngout (bool)
                Optional, write a file to disk?
            fname (str)
                Output filename, required only if pngout=True            
        """

        # The simple figsize format is OK for static images:
        # fig = plt.figure(figsize=(20,9))

        # If saved files are to be used for animations with ffmpeg, row
        # count must be an even number. This hack ensures that.
        fig = plt.figure()
        DPI = fig.get_dpi() # dots per inch of your display
        fig.set_size_inches(1200.0/float(DPI),610.0/float(DPI))

        subplots = (121, 122)

        # set up the left subplot: Rotated Disk - FACE ON
        ax0 = plt.subplot(121)

        # plot the particle density
        # can modify bin number (bin =100 smoothest)
        h0 = ax0.hist2d(rn[0,:], rn[1,:], bins=nbins, norm=LogNorm(), cmap='magma')
        fig.colorbar(h0[3], ax=ax0)

        # Add axis labels
        fontsize = 18
        ax0.set_xlabel('x (kpc)', fontsize=fontsize)
        ax0.set_ylabel('y (kpc)', fontsize=fontsize)

        #set axis limits
        if lim is None:
            ax0.axis('equal')
        else:
            ax0.set_ylim(-lim,lim)
            ax0.set_xlim(-lim,lim)

        # make the contour plot
        # x pos, y pos, contour res, contour res, axis, colors for contours.
        level_vals = [0.68, 0.8, 0.9, 0.95, 0.99]
        colors = ['red','orange', 'yellow', 'orange', 'yellow']
        self.density_contour(rn[0,:], rn[1,:], 80, 80, level_vals, ax=ax0, colors=colors)

        # set up the right subplot: Rotated Disk - EDGE ON
        ax1 = plt.subplot(122)

        # plot the particle density
        # can modify bin number (bin =100 smoothest)
        h1 = ax1.hist2d(rn[0,:], rn[2,:], bins=nbins, norm=LogNorm(), cmap='magma')
        fig.colorbar(h1[3], ax=ax1)

        # Add axis labels
        ax1.set_xlabel('x (kpc)', fontsize=fontsize)
        ax1.set_ylabel('z (kpc)', fontsize=fontsize)

        #set axis limits
        if lim is None:
            ax1.axis('equal')
        else:
            ax1.set_ylim(-lim,lim)
            ax1.set_xlim(-lim,lim)

        self.density_contour(rn[0,:], rn[2,:], 80, 80, level_vals, ax=ax1, colors=colors)

        # adjusting tick label font size led to all sorts of problems
        # so commented out for now

        # label_size = 18
        # matplotlib.rcParams['xtick.labelsize'] = label_size 
        # matplotlib.rcParams['ytick.labelsize'] = label_size

        fig.suptitle(f'{galname}, t = {t:5.2f} Gyr', fontsize=26, weight='bold')

        # Save to a file
        if pngout:
            plt.savefig(fname, dpi='figure');   
            
    def plot_phase(self, rn, vn, R, Vcirc, galname, t, xlim=20, ylim=200, 
                    bins=500, colorbar=True, pngout=False, fname=None):
        """
        """

        # Make a phase diagram
        # Disk Velocity Field edge on.

        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(1000.0/float(DPI),700.0/float(DPI))

        ax = plt.subplot(111)

        # looking at galaxy edge on along x axis, vy is line of sight velocity

        plt.hist2d(rn[0,:], vn[1,:], bins=bins, norm=LogNorm())
        if colorbar:
            plt.colorbar()

        # Add the circular velocity
        plt.plot(R, Vcirc, color="red")
        plt.plot(-R, -Vcirc, color="red")

        # Add axis labels
        plt.xlabel('x (kpc)', fontsize=22)
        plt.ylabel('Velocity Y Direction (km/s)', fontsize=22)

        plt.xlim(-xlim, xlim)
        plt.ylim(-ylim, ylim)

        plt.title(f'{galname}, t = {t:5.2f} Gyr', fontsize=26, weight='bold')  

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

    def mass_profiles(self, radii, mass_profiles, t, figsize=(18,8), ylim=(1e7, 1e12), 
                    pngout=False, fname=None):
        """
        The `MassProfile` class has methods `mass_enclosed()` and `mass_enclosed_total()` 
        which accept arrays of radii. This gives us most of what we need for plotting.

        Args:

        """

        assert len(mass_profiles) == 3 # MW, M31, M33 - matching the 3 subplots
        fig = plt.figure(figsize=figsize)
        subplots = (131, 132, 133)

        for i in range(len(mass_profiles)):
            # set up this subplot
            ax = plt.subplot(subplots[i])
            gname, mp = mass_profiles[i]
            
            # add the curves
            ax.semilogy(radii, mp.mass_enclosed_total(radii), 'b-', lw=3, 
                        label='Total')
            ax.semilogy(radii, mp.mass_enclosed(radii, 1), 'r:', lw=3, 
                        label='Halo particles')
            ax.semilogy(radii, mp.mass_enclosed(radii, 2), 'g--', lw=3, 
                        label='Disk particles')
            ax.semilogy(radii, mp.mass_enclosed(radii, 3), 'm-.', lw=3, 
                        label='Bulge particles')

            #adjust tick label font size
            label_size = 16
            rcParams['xtick.labelsize'] = label_size 
            rcParams['ytick.labelsize'] = label_size

            # Add labels and subplot title
            ax.set_xlabel('Radius (kpc)', fontsize=20)
            if i == 0: # left subplot only
                ax.set_ylabel(r'Mass enclosed ($M_\odot$)', fontsize=22)
                ax.legend(loc='lower right',fontsize='xx-large', shadow=True)
            ax.set_title(gname[:-4], fontsize=24)

            #set axis limits
            ax.set_ylim(ylim[0], ylim[1])

        # Overall title
        fig.suptitle(f'Mass Profiles by Particle Type, t={t:.2f} Gyr', y=1.0, 
                    fontsize=24, weight='bold')

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

    def rotation_curves(self, radii, mass_profiles, fitted_a, t, 
                 figsize=(18,8), ylim=(1e7, 1e12), pngout=False, fname=None):
        """
        """
        
        assert len(mass_profiles) == 3
        fig = plt.figure(figsize=figsize)
        subplots = (131, 132, 133)

        for i in range(len(mass_profiles)):
            ax = plt.subplot(subplots[i])
            gname, mp = mass_profiles[i]
            
            # the a value is from scipy.optimize, above
            a_opt = fitted_a[gname]
            
            # add the curves
            ax.plot(radii, mp.circular_velocity_total(radii), 'b-', lw=3, 
                        label='Total')
            ax.plot(radii, mp.circular_velocity(radii, 1), 'r:', lw=3, 
                        label='Halo particles')
            ax.plot(radii, mp.circular_velocity(radii, 2), 'g--', lw=3, 
                        label='Disk particles')
            
            # M33 has no bulge, but include it just to get a legend
            # The other subplots have no empty areas
            ax.plot(radii, mp.circular_velocity(radii, 3), 'm-.', lw=3, 
                        label='Bulge particles')
            
            # add the Hernquist profile, use our best-fit a 
            ax.plot(radii, mp.circular_velocity_hernquist(radii, a_opt), 'c-', lw=2, 
                        label='Hernquist profile')

            #adjust tick label font size
            label_size = 16
            rcParams['xtick.labelsize'] = label_size 
            rcParams['ytick.labelsize'] = label_size

            # Add labels
            ax.set_xlabel('Radius (kpc)', fontsize=20)
            if i == 0: # left subplot only
                ax.set_ylabel(r'Circular velocity (km/s)', fontsize=22)
            if i == 2: # right subplot only
                ax.legend(fontsize='xx-large', shadow=True)
            ax.set_title(gname[:-4], fontsize=24)

            #set axis limits
            ax.set_ylim(0, 270)

        fig.suptitle(f'Rotation curves by Particle Type, t={t:.2f}', y=1.0, 
                    fontsize=24, weight='bold')

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

    def plot_h_r(self, radii, h, xmax=350, figsize=(8,6), ylabel='$h$ within radius',
                pngout=False, fname=None):
        fontsize = 24
        fig = plt.figure(figsize=figsize)
        plt.plot(radii, norm(h, axis=1), 'b-', lw=3)
        plt.xlabel('r (kpc)', fontsize=fontsize)
        plt.ylabel(ylabel, fontsize=fontsize)
        plt.xlim(0, xmax),
        # plt.ylim(0, 3000)

        label_size = 16
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size

        plt.tight_layout()
        if pngout:
            plt.savefig(fname, rasterized=True, dpi=350);

    def plot_theta_phi(self, radii, theta, phi, xmax=350, figsize=(8,6),
                pngout=False, fname=None):
        fontsize = 24
        fig = plt.figure(figsize=figsize)
        plt.plot(radii, theta, 'r:', lw=3, label=r'$\theta$ (polar)')
        plt.plot(radii, phi, 'b-', lw=3, label=r'$\phi$ (azimuthal)')
        plt.xlabel('r (kpc)', fontsize=fontsize)
        plt.ylabel('$\hat{L}$ angles (deg)', fontsize=fontsize)
        plt.xlim(0, xmax)
        # plt.ylim(40, 90)

        label_size = 16
        matplotlib.rcParams['xtick.labelsize'] = label_size 
        matplotlib.rcParams['ytick.labelsize'] = label_size

        plt.legend(fontsize='xx-large', shadow=True)
        plt.tight_layout()
        if pngout:
            plt.savefig(fname, rasterized=True, dpi=350);

    def plot_v_sigma(self, xbins, means, sigmas, xlim=None, ylim1=None, ylim2=None, 
                xlabel=None, title=None, figsize=(8,6), pngout=False, fname=None):
        
        if xlabel is None:
            xlabel = 'y (kpc)'

        fig, ax1 = plt.subplots(figsize=figsize)

        color = 'red'
        ax1.plot(xbins, means, 'r.', label='radial velocity')

        fontsize = 24
        ax1.set_xlabel(xlabel, fontsize=fontsize)
        ax1.set_ylabel('Mean radial velocity (km/s)', color=color, fontsize=fontsize)
        ax1.tick_params(axis='y', labelcolor=color)
        if xlim is not None:
            ax1.set_xlim(xlim[0], xlim[1])
        if ylim1 is not None:
            ax1.set_ylim(ylim1[0], ylim1[1])

        ax2 = ax1.twinx()  # shares the same x-axis as ax1

        color = 'blue'
        ax2.plot(xbins, sigmas, 'b+', label='velocity dispersion')
        ax2.set_ylabel('Velocity dispersion (km/s)', color=color, fontsize=fontsize)
        ax2.tick_params(axis='y', labelcolor=color)
        if ylim2 is not None:
            ax2.set_ylim(ylim2[0], ylim2[1])

        # plt.legend()

        if title is not None:
            ax1.set_title(title, fontsize=24)

        #adjust tick label font size
        label_size = 18
        rcParams['xtick.labelsize'] = label_size 
        rcParams['ytick.labelsize'] = label_size

        fig.tight_layout()  
        if pngout:
            plt.savefig(fname, rasterized=True, dpi=350);
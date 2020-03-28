# import modules
import numpy as np
# from numpy.linalg import norm
import scipy.optimize as so

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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
            level_vals : list of int
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
            
    def plot_phase(self, rn, vn, R, Vcirc, galname, t, xlim=20, ylim=200, pngout=False, fname=None):
        """
        """

        # Make a phase diagram
        # Disk Velocity Field edge on.

        fig = plt.figure()
        DPI = fig.get_dpi()
        fig.set_size_inches(1000.0/float(DPI),700.0/float(DPI))

        ax = plt.subplot(111)

        # looking at galaxy edge on along x axis, vy is line of sight velocity

        plt.hist2d(rn[0,:], vn[1,:], bins=150, norm=LogNorm())
        plt.colorbar()

        # Add the circular velocity
        plt.plot(R, Vcirc, color="red")
        plt.plot(-R, -Vcirc, color="red")

        # Add axis labels
        plt.xlabel('x (kpc)', fontsize=22)
        plt.ylabel('Velocity Y Direction (km/s)', fontsize=22)

        plt.xlim(-xlim, xlim)
        plt.ylim(-ylim, ylim)

        # Save file
        if pngout:
            plt.savefig(fname, dpi='figure');   

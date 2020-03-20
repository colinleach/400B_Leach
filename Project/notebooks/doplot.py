# import modules
import numpy as np
from numpy.linalg import norm
import astropy.units as u
from astropy.constants import G
from pathlib import Path

# import plotting modules
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

# my modules
from galaxy.galaxy import Galaxy
from galaxy.centerofmass import CenterOfMass
from galaxy.massprofile import MassProfile
from galaxy.timecourse import TimeCourse
from galaxy.plots import Plots

p = Plots()

gname = 'MW'

for snap in np.arange(0, 300):
    fname = f'{gname}_density_{snap:03}.png'
    print(fname)
    if not Path(fname).is_file():
        print(snap)
        try:
            gal = Galaxy(gname, snap, usesql=True, ptype=2)
            t = gal.time.value / 1000
        except TypeError:
            gal = Galaxy(gname, snap, datadir='/home/colin/zcode/UAclasses/400b/HighRes', ptype=2)
            t = gal.time.value / 1000

        com = CenterOfMass(gal)

        tc = TimeCourse(usesql=True)
        com_xyz, com_vxyz = tc.get_one_com(gname, snap)

        gal_xyzD, gal_vxyzD = com.center_com(com_xyz, com_vxyz)

        # determine the rotated velocity vectors
        rn, vn = com.rotate_frame(com_p=com_xyz, com_v=com_vxyz)

        p.plot_density(rn, gname, t, pngout=True, snap=snap, lim=60)
        plt.close('all')
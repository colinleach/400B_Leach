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

def make_plot(gname, snap, lim):
    try:
        gal = Galaxy(gname, snap, usesql=True, ptype=2)
        t = gal.time.value / 1000
    except TypeError:
        gal = Galaxy(gname, snap, datadir=datadir, ptype=2)
        t = gal.time.value / 1000

    com = CenterOfMass(gal)

    tc = TimeCourse(usesql=True)
    com_xyz, com_vxyz = tc.get_one_com(gname, snap)

    gal_xyzD, gal_vxyzD = com.center_com(com_xyz, com_vxyz)

    # determine the rotated velocity vectors
    rn, vn = com.rotate_frame(com_p=com_xyz, com_v=com_vxyz)

    p.plot_density(rn, gname, t, pngout=True, snap=snap, lim=60)
    plt.close('all')

p = Plots()

limits = {'MW': (50, 80),
        'M31': (50, 80), 
        'M33': (30, 100)}

cmd = ''

datadir = Path.home() / 'HighRes'

for gname in ('MW', 'M31', 'M33'):
    print(gname)
    for snap in np.arange(0, 300):
        print(snap, end=' ')
        lim = limits[galname][0]
        make_plot(gname, snap, lim, group='early')
    cmd += f'ffmpeg -r 10  -start_number 0  -s 1920x1080 -i {gname}_density_early_%03d.png'
    cmd += f' -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p {gname}_early.mp4\n'

    for snap in np.arange(290, 802):
        print(snap, end=' ')
        lim = limits[galname][1]
        make_plot(gname, snap, lim, group='late')
    cmd += f'ffmpeg -r 10  -start_number 0  -s 1920x1080 -i {gname}_density_late_%03d.png'
    cmd += f' -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p {gname}_late.mp4\n'


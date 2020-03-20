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

def make_plot(gname, snap, lim, fname):
    try:
        gal = Galaxy(gname, snap, usesql=True, ptype=2)
        t = gal.time.value / 1000
    except TypeError:
        gal = Galaxy(gname, snap, datadir=datadir, ptype=2)
        t = gal.time.value / 1000

    com = CenterOfMass(gal)

    tc = TimeCourse(usesql=True)
    com_xyz, com_vxyz = tc.get_one_com(gname, snap)

    # gal_xyzD, gal_vxyzD = com.center_com(com_xyz, com_vxyz)

    # determine the rotated velocity vectors
    rn, _ = com.rotate_frame(com_p=com_xyz, com_v=com_vxyz)

    p.plot_density(rn, gname, snap, t, pngout=True, lim=60, fname=fname)
    plt.close('all')

p = Plots()

limits = {'MW': (50, 80),
        'M31': (50, 80), 
        'M33': (30, 100)}

cmd = ''

datadir = Path.home() / 'HighRes'
cmdfile = 'make_densities.sh'

with open(cmdfile, 'w') as fp:
    fp.write(cmd)

for gname in ('MW', 'M31', 'M33'):
    print(gname)
    group = 'early'
    for snap in np.arange(0, 300):
        print(snap, end=' ')
        lim = limits[gname][0]
        fname = f'png_files/{gname}_density_{group}_{snap:03}.png'
        make_plot(gname, snap, lim, fname=fname)
    cmd += f'ffmpeg -r 10  -start_number 290  -s 1920x1080'
    cmd += f' -i png_files/{gname}_density_early_%03d.png'
    cmd += f' -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p {gname}_early.mp4\n'
    with open(cmdfile, 'w') as fp:
        fp.write(cmd)

    for snap in np.arange(290, 802):
        print(snap, end=' ')
        group = 'late'
        lim = limits[gname][1]
        fname = f'png_files/{gname}_density_{group}_{snap:03}.png'
        make_plot(gname, snap, lim, fname=fname)
    cmd += f'ffmpeg -r 10  -start_number 290  -s 1920x1080'
    cmd += f' -i png_files/{gname}_density_late_%03d.png'
    cmd += f' -vcodec libx264 -vf fps=25 -crf 25  -pix_fmt yuv420p {gname}_late.mp4\n'
    with open(cmdfile, 'w') as fp:
        fp.write(cmd)


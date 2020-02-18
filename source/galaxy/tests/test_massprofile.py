import numpy as np
# import pandas as pd
from numpy.linalg import norm
import astropy.units as u

import pytest
from pytest import approx

from galaxy import galaxies
from galaxy import massprofile

g = galaxies.Galaxies()

def test_mp_create():
    mp = massprofile.MassProfile(g.galaxies['MW_000'])
    expected = np.array([-2.07, 2.95, -1.45])
    assert mp.com_p.value == approx(expected)

mp = massprofile.MassProfile(g.galaxies['MW_000'])

def test_mass_enclosed():
    mass = mp.mass_enclosed([1, 5, 20]*u.kpc, ptype=2)
    expected = [1.1206e+10, 3.8671e+10, 7.1190e+10]
    assert mass.value == approx(expected)

def test_mass_enclosed_total():
    mass = mp.mass_enclosed([1, 5, 20]*u.kpc)
    expected = [1.57065155e+10, 5.60813060e+10, 1.99080195e+11]
    assert mass.value == approx(expected)
   
def test_halo_mass():
    mass = mp.halo_mass()
    expected = 1.974925e12
    assert mass.value == approx(expected)
   
def test_hernquist_mass():
    mass = mp.hernquist_mass([1, 5, 20]*u.kpc, 60*u.kpc)
    expected = [5.30751142e+08, 1.16859467e+10, 1.23432812e+11]
    assert mass.value == approx(expected)
   

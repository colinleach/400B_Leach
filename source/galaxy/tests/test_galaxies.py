import numpy as np
import pandas as pd

import pytest
from pytest import approx

from galaxy import galaxies

def test_galaxies_create():
    g = galaxies.Galaxies()
    assert g.filenames == ['MW_000', 'M31_000', 'M33_000']

g = galaxies.Galaxies()

@pytest.mark.parametrize('gal, size', [('MW_000', 135000),
                                        ('M31_000', 189000),
                                        ('M33_000', 14300),])
def test_galaxies_data(gal, size):
    assert g.galaxies[gal].data.shape[0] == size

def test_count_pivot():
    cp = g.get_counts_pivot()
    assert cp.shape == (4, 4)
    assert cp['1 Halo']['All'] == 105000

def test_masses_pivot():
    mp = g.get_masses_pivot()
    assert mp.shape == (4, 4)
    assert mp['1 Halo']['All'] == 4082418000000

def test_get_coms():
    coms = g.get_coms()
    colnames = ['name', 'ptype', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'R', 'V']
    assert coms.colnames == colnames
    V_vals = np.array([14.14,6.53,10.2,115.61,113.63,110.05,181.72,180.45])
    assert coms['V'].value[:-1] == approx(V_vals)

def test_separations():
    seps = g.separations('MW_000', 'M31_000')
    assert seps['r'].value == approx(769.1)
    assert seps['vel_mag'].value == approx(117.74)
    assert seps['v_radial'].value == approx(-115.85)
    assert seps['v_tan_mag'].value == approx(21.36)
    
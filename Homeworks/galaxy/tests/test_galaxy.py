import numpy as np
import pandas as pd

from pytest import approx

from galaxy import galaxy


def test_galaxy_init():
    MW = galaxy.Galaxy('MW', snap=0)
    assert MW.filename == 'MW_000.txt'

MW = galaxy.Galaxy('MW', snap=0)

def test_MW_data():
    assert MW.data.shape == (135000, )

def test_filter_by_type():
    type = 2 # just disk particles
    disk_particles = MW.filter_by_type(2)
    assert disk_particles.shape == (75000,)

def test_single_particle_properties():
    type = 2
    num = 99
    pos, v, m = MW.single_particle_properties(particle_num=num, type=type)

    # get the corresponding particle from the full dataset
    subset = MW.filter_by_type(type)
    raw_prop = subset[num]

    # compare 2 methods of calculating
    raw_pos = np.sqrt(raw_prop['x']**2 + raw_prop['y']**2 + raw_prop['z']**2)
    assert pos.value == approx(raw_pos, rel=1e-4)
    raw_v = np.sqrt(raw_prop['vx']**2 + raw_prop['vy']**2 + raw_prop['vz']**2)
    assert v.value == approx(raw_v, rel=1e-4)

    assert m.value == approx(raw_prop['m'] * 1e10)

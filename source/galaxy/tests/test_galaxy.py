import numpy as np
import pandas as pd

import pytest
from pytest import approx

from galaxy import galaxy


@pytest.mark.parametrize('gal, fname', [('MW', 'MW_000.txt'),
                                        ('M31', 'M31_000.txt'),
                                        ('M33', 'M33_000.txt'),])
def test_galaxy_init(gal, fname):
    gal_obj = galaxy.Galaxy(gal, snap=0)
    assert gal_obj.filename == fname

galaxies = {
    'MW': galaxy.Galaxy('MW', snap=0),
    'M31': galaxy.Galaxy('M31', snap=0),
    'M33': galaxy.Galaxy('M33', snap=0),
}

@pytest.mark.parametrize('gal, size', [('MW', 135000),
                                        ('M31', 189000),
                                        ('M33', 14300),])
def test_MW_data(gal, size):
    assert galaxies[gal].data.shape[0] == size

@pytest.mark.parametrize('gal, size', [('MW', (50000, 75000, 10000)),
                                        ('M31', (50000, 120000, 19000)),
                                        ('M33', (5000, 9300, 0)),])
def test_filter_by_type(gal, size):
    for ptype in (1,2,3):
        particles = galaxies[gal].filter_by_type(ptype)
        assert particles.shape[0] == size[ptype - 1]

@pytest.mark.parametrize('gal', ('MW', 'M31', 'M33'))
def test_single_particle_properties(gal):
    for ptype in (None,1,2,3):
        num = 99
        try:
            pos, v, m = galaxies[gal].single_particle_properties(particle_num=num, 
                                                                particle_type=ptype)
        except IndexError: # no bulge particles in M33
            continue # skip the rest for this galaxy/ptype combination

        # get the corresponding particle from the full dataset
        if ptype is None:
            subset = galaxies[gal].data
        else:
            subset = galaxies[gal].filter_by_type(ptype)
        raw_prop = subset[num]

        # compare 2 methods of calculating
        raw_pos = np.sqrt(raw_prop['x']**2 + raw_prop['y']**2 + raw_prop['z']**2)
        assert pos.value == np.around(raw_pos, 3)
        raw_v = np.sqrt(raw_prop['vx']**2 + raw_prop['vy']**2 + raw_prop['vz']**2)
        assert v.value == np.around(raw_v, 3)
        assert m.value == approx(raw_prop['m'] * 1e10)

@pytest.mark.parametrize('gal, size', [('MW', (135000, 50000, 75000, 10000)),
                                        ('M31', (189000, 50000, 120000, 19000)),
                                        ('M33', (14300, 5000, 9300, 0)),])
def test_all_particle_properties(gal, size):
    for ptype in (None,1,2,3):
        qtprops = galaxies[gal].all_particle_properties(particle_type=ptype)
        assert qtprops.colnames == ['type', 'm', 'pos', 'v']
        if ptype is None:
            assert len(qtprops) == size[0]
        else:
            assert len(qtprops) == size[ptype]

@pytest.mark.parametrize('gal, counts', [('MW', [50000, 75000, 10000]),
                                        ('M31', [50000, 120000, 19000]),
                                        ('M33', [5000, 9300, 0]),])
def test_all_component_counts(gal, counts):
    assert galaxies[gal].all_component_counts() == counts
         
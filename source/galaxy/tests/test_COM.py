import numpy as np
import pandas as pd
from numpy.linalg import norm

import pytest
from pytest import approx

from galaxy import galaxies
from galaxy import centerofmass

g = galaxies.Galaxies()

def test_com_create():
    com = centerofmass.CenterOfMass(g.galaxies['MW_000'], 2)
    assert com.xyz.shape == (3, 75000)

MW_2 = centerofmass.CenterOfMass(g.galaxies['MW_000'], 2)

def test_com_define():
    com_xyz = MW_2.com_define(MW_2.xyz, MW_2.m)
    expected = np.array([-0.93398577,  2.40892628, -1.42421335])
    assert com_xyz == approx(expected)

def test_com_p():
    com_p = MW_2.com_p()
    expected = [-2.07, 2.95, -1.45]
    assert com_p.value == approx(expected)

def test_com_v():
    com_p = MW_2.com_p()
    com_v = MW_2.com_v(com_p)
    expected = [0.94, 6.32, -1.35]
    assert com_v.value == approx(expected)
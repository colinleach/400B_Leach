# ParticleProperties.py
#
# Colin Leach
#
# All the interesting code is in the Galaxy class (../galaxy/galaxy.py)
# This file is just syntactic sugar to comply with homework requirements

# standard imports
import numpy as np
import astropy.units as u

# my module import currently relies on symlink (ugly):
#   `ln -s ../galaxy/galaxy.py galaxy.py`
from ReadFile import Read, ReadGalaxy


def ParticleInfo(galaxy_name, particle_type, particle_number):
    "Wrapper around the Galaxy.single_particle_properties() method"

    gal = ReadGalaxy(galaxy_name)
    return gal.single_particle_properties(particle_type, particle_number)


if __name__ == '__main__':
    ptypes = {'DM': 1, 'Disk': 2, 'Bulge': 3}

    # we want the 100th disk particle in the Milky Way:
    galaxy_name = 'MW'
    ptype_name = 'Disk'
    ptype = ptypes[ptype_name]
    pnum = 99  # zero-based index

    pos, v, m = ParticleInfo(galaxy_name,
                             particle_type=ptype,
                             particle_number=pnum)

    print(f'{ptype_name} particle # {pnum + 1} in {galaxy_name}:')
    print(f'Distance from CoM: {pos},  velocity: {v},  mass: {m}')
    print(f'In light years, distance is {np.around(pos.to(u.lyr), 3)}')

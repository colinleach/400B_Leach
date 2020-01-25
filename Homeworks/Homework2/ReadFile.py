# ReadFile.py
#
# Colin Leach
#
# All the interesting code is in the Galaxy class (../galaxy/galaxy.py)
# This file is just syntactic sugar to comply with homework requirements

# import currently relies on symlink: `ln -s ../galaxy/galaxy.py galaxy.py`
from galaxy import Galaxy


def Read(galaxy_name='MW', snap=0):
    "Returns the 3-tuple asked for in Hw2"

    gal = Galaxy(galaxy_name, snap)
    return gal.time, gal.particle_count, gal.data


def ReadGalaxy(galaxy_name='MW', snap=0):
    "Returns the Galaxy object, with various useful methods"

    return Galaxy(galaxy_name, snap)


if __name__ == '__main__':
    time, count, data = Read()
    print(f'time: {time},  particle count: {count},  data shape: {data.shape}')

# standard Python imports
from pathlib import Path

# scientific package imports
import numpy as np
from numpy.linalg import norm

import astropy.units as u
from astropy.table import QTable

import pandas as pd


class Galaxy():
    """
    Args: 
        name (str): 
            short name used in filename of type 'name_000.txt', eg 'MW', 'M31'.

    Kwargs:
        snap (int): 
            Snap number, equivalent to time elapsed. Zero is starting conditions.
        datadir (str):
            Directory to search first for the required file. Optional, and a 
            default list of locations will be searched.

    Class attributes:
        path (`pathlib.Path` object): 
            directory containing the data file
        filename (str):
            in `name_snap.txt` format, something like 'MW_000.txt'
        data (np.ndarray):
            type, mass, position_xyz, velocity_xyz for each particle

    A class to find, read and manipulate files for a single galaxy.
    """

    def __init__(self, name, snap=0, datadir=None):
        "Initial setup. Currently it calls read_file(), but this may change."

        self.name = name
        self.snap = snap

        # We can probably make some assumptions about the filename:
        self.filename = f"{name}_{snap:03}.txt"
        self.filepath = self.get_filepath(datadir)
        self.read_file()

    def get_filepath(self, datadir):
        """
        Args:
            datadir (str): path to search first for the required file

        Returns:   
            `pathlib.Path` object. A directory containing the file.

        Raises:
            FileNotFoundError

        Pretty boring housekeeping code, but may make things more resilient.
        """

        # helper function, returns (file found at this location?) True/False
        def has_file(dirname):
            if dirname is not None and \
                    Path(dirname).exists() and \
                    Path(dirname).is_dir():
                if (Path(dirname) / self.filename).exists():
                    return True
            return False

        if has_file(datadir):
            return Path(datadir)  # whoever called this knew the layout

        # Now we aren't sure where the file is so may need to search for it
        # it may be different on different systems, e.g. nimoy vs laptop
        # try the search path:
        #   [datadir, '.', '..', '../..', '~', '~/galaxydata']
        pwd = Path.cwd()
        home = Path.home()
        searchpath = [pwd, ]
        searchpath += [pwd.parents[n] for n in range(4)]
        searchpath += [home, home / 'galaxydata']

        for p in searchpath:
            if has_file(p):
                return p  # happy result, we got a valid path

        # Raise an error if not found
        pathstrings = f'{datadir}, ' if datadir is not None else ''
        pathstrings += ', '.join([sp.as_posix() for sp in searchpath])
        msg = f'{self.filename}: Not found in {pathstrings}'
        raise FileNotFoundError(msg)

    def read_file(self):
        """
        Read in a datafile in np.ndarray format, store in `self.data`.

        Requires:
            `self.path` and `self.filename` are already set.

        Returns: 
            nothing
        """

        fullname = self.filepath / self.filename

        # get header data
        with open(fullname) as file:
            _, value = file.readline().split()
            self.time = float(value) * 10.0 * u.Myr
            label, value = file.readline().split()
            self.particle_count = int(value)

        # get the big array of values
        self.data = np.genfromtxt(
            fullname, dtype=None, names=True, skip_header=3)

    def filter_by_type(self, type, dataset=None):
        """
        Args:
            type (int): for particles, 1=DM, 2=disk, 3=bulge 
    
        Kwargs:
            dataset (np.ndarray): a starting dataset other than self.data

        Returns: 
            np.ndarray: subset data
        """

        if dataset is None:
            dataset = self.data

        index = np.where(dataset['type'] == type)
        return dataset[index]

    def single_particle_properties(self, type=None, particle_num=0):
        """
        Kwargs:
            type (int): 
                a subset of the data filtered by 1=DM, 2=disk, 3=bulge
            particle_num (int): 
                zero-based index to an array of particles

        returns: 
            3-tuple of
                Euclidean distance from CoM (kpc),
                Euclidean velocity magnitude (km/s),
                particle mass (M_sun)
        """

        # The next bit will throw IndexError if particle_num invalid
        # Be ready to catch this
        if type is None:  # all types accepted
            particle = self.data[particle_num]
        else:
            particle = self.filter_by_type(type)[particle_num]

        # mass:
        m = particle['m'] * 1e10 * u.Msun

        # Euclidean distance from galactic CoM:
        pos_xyz = np.array([particle['x'], particle['y'], particle['z']])
        pos_mag = norm(pos_xyz) * u.kpc

        # velocity:
        v_xyz = np.array([particle['vx'], particle['vy'], particle['vz']])
        v_mag = norm(v_xyz) * u.km / u.s

        return np.around(pos_mag, 3), np.around(v_mag, 3), m

    def all_particle_properties(self, type=None):
        """
        Kwargs:
            type (int): 
                a subset of the data filtered by 1=DM, 2=disk, 3=bulge

        Returns:
            QTable: The full list with units, optionally filtered by type.
        """

        if type is None:  # all types accepted
            dataset = self.data
        else:
            dataset = self.filter_by_type(type)

        # mass:
        m = dataset['m'] * 1e10 * u.Msun

        # Pythagorean distance from galactic CoM:
        pos_xyz = np.array([dataset['x'], dataset['y'], dataset['z']])
        pos_mag = norm(pos_xyz, axis=0) * u.kpc

        # velocity:
        v_xyz = np.array([dataset['vx'], dataset['vy'], dataset['vz']])
        v_mag = norm(v_xyz, axis=0) * u.km / u.s

        # construct and return a QTable with units
        t = QTable()
        t['type'] = dataset['type']
        t['m'] = np.around(m)
        t['pos'] = np.around(pos_mag, 3)
        t['v'] = np.around(v_mag, 3)
        return t

    # ________________________________________________________________
    #
    # define some getters which may turn out to be useful, perhaps

    def get_array(self):
        """
        Returns: 
            data in `np.ndarray` format

        Pretty superfluous in Python (which has no private class members)
        """

        return self.data

    def get_df(self):
        """Returns:
            data as pandas dataframe
        """

        return pd.DataFrame(self.data)

    def get_qtable(self):
        """Returns:
            data as astropy QTable, with units
        """

        t = QTable(self.data)

        # add appropriate units
        t['m'] = np.around(t['m'] * 1e10) * u.Msun
        for col in ('x', 'y', 'z'):
            t[col] *= u.kpc
        for col in ('vx', 'vy', 'vz'):
            t[col] *= (u.km / u.s)

        return t
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
    A class to find, read and manipulate files for a single galaxy.

    Args:
        name (str):
            short name used in filename of type 'name_000.txt', eg 'MW', 'M31'.
        snap (int):
            Snap number, equivalent to time elapsed. Zero is starting conditions.
        datadir (str):
            Directory to search first for the required file. Optional, and a
            default list of locations will be searched.
        usesql (bool):
            If True, data will be taken from a PostgreSQL database instead of
            text files.
        ptype (int):
            Optional. Restrict data to this particle type, for speed. 
            Only valid with usesql=True.
        stride (int):
            Optional. For stride=n, get every nth row in the table.
            Only valid with usesql=True.

    Class attributes:
        filepath (`pathlib.Path` object):
            directory containing the data file
        filename (str):
            in `name_snap.txt` format, something like 'MW_000.txt'
        data (np.ndarray):
            type, mass, position_xyz, velocity_xyz for each particle
    """

    def __init__(self, name, snap=0, datadir=None, usesql=False, ptype=None, stride=1):
        "Initial setup. Currently it calls read_file(), but this may change."

        self.name = name
        self.snap = snap

        if usesql:
            self.read_db(ptype, stride)
        else:
            # We can probably make some assumptions about the filename:
            self.filename = f"{name}_{snap:03}.txt"
            self.filepath = self.get_filepath(datadir)
            self.read_file()

    def read_db(self, ptype, stride):
        """
        Get relevant data from a PostgreSQL database and format it to be 
        identical to that read from test files.

        Args:
            ptype (int):
                Optional. Restrict data to this particle type.
            stride (int):
                Optional. For stride=n, get every nth row in the table.

        Changes:
            `self.time`, `self.particle_count` and `self.data` are set.

        Returns: nothing
        """

        from galaxy.db import DB

        db = DB()
        cur = db.get_cursor()

        # set the elapsed time
        sql_t = f"SELECT time FROM simdata WHERE galname='{self.name}'"
        sql_t += f" and snap={self.snap} LIMIT 1"
        cur.execute(sql_t)
        time = cur.fetchone()
        self.time = time[0]

        # set the bulk of the data
        colheads = ','.join(['type','m','x','y','z','vx','vy','vz'])
        if stride > 1:
            sql_d = f"SELECT {colheads}, ROW_NUMBER() OVER () as rn from simdata"
        else:
            sql_d = f"SELECT {colheads} from simdata"
        sql_d += f"  where galname='{self.name}' and snap={self.snap}"
        if ptype:
            sql_d += f" and type={ptype}"
        if stride > 1:
            sql_d = f"SELECT {colheads} from ( {sql_d} ) as t where rn % {stride} = 0" 

        dtype=[('type', 'uint8'), ('m', '<f4'), ('x', '<f4'), ('y', '<f4'), ('z', '<f4'), 
                ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        cur.execute(sql_d)
        self.data = np.array(cur.fetchall(), dtype=dtype)
        self.particle_count = len(self.data)

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

        Changes:
            `self.time`, `self.particle_count` and `self.data` are set.

        Returns: nothing
        """

        fullname = self.filepath / self.filename

        # get header data
        with open(fullname) as file:
            _, value = file.readline().split()
            self.time = float(value) * u.Myr  # corrected previous unit error
            label, value = file.readline().split()
            self.particle_count = int(value)

        # get the big array of values
        self.data = np.genfromtxt(
            fullname, dtype=None, names=True, skip_header=3)

    # --------------------------------------------------------------------
    # A couple of small utility functions to interconvert particle type
    # --------------------------------------------------------------------

    def type2name(self, particle_type):
        """
        Args: particle_type (int): valid values are 1, 2, or 3

        Returns: typename (str): 'DM', 'disk' or 'bulge'
        """

        typenames = {1: 'DM', 2: 'disk', 3: 'bulge'}
        return typenames[particle_type]

    def name2type(self, typename):
        """
        Args: typename (str): valid values are 'DM', 'disk' or 'bulge'

        Returns: particle_type (int): 1, 2, or 3 as used in data files
        """

        types = {'DM': 1, 'disk': 2, 'bulge': 3}
        return types[typename]

    # --------------------------------------------------------------------
    # Subset the data by type
    # --------------------------------------------------------------------

    def filter_by_type(self, particle_type, dataset=None):
        """
        Subsets the data to a single particle type.

        Args:
            particle_type (int): for particles, 1=DM, 2=disk, 3=bulge
            dataset (array including a type column): defaults to self.data

        Kwargs:
            dataset (np.ndarray): optionally, a starting dataset other than self.data

        Returns: np.ndarray: subset data
        """

        if dataset is None:
            dataset = self.data

        index = np.where(dataset['type'] == particle_type)
        return dataset[index]

    # --------------------------------------------------------------------
    # Methods to calculate properties of an existing dataset
    # --------------------------------------------------------------------

    def single_particle_properties(self, particle_type=None, particle_num=0):
        """
        Calculates distance from the origin and magnitude of the velocity.

        Kwargs:
            particle_type (int):
                a subset of the data filtered by 1=DM, 2=disk, 3=bulge
            particle_num (int):
                zero-based index to an array of particles

        returns:
            3-tuple of
                Euclidean distance from origin (kpc),
                Euclidean velocity magnitude (km/s),
                particle mass (M_sun)
        """

        # The next bit will throw IndexError if particle_num invalid
        # Be ready to catch this
        if particle_type is None:  # all types accepted
            particle = self.data[particle_num]
        else:
            particle = self.filter_by_type(particle_type)[particle_num]

        # mass:
        m = particle['m'] * 1e10 * u.Msun

        # Euclidean distance from origin:
        pos_xyz = np.array([particle['x'], particle['y'], particle['z']])
        pos_mag = norm(pos_xyz) * u.kpc

        # velocity:
        v_xyz = np.array([particle['vx'], particle['vy'], particle['vz']])
        v_mag = norm(v_xyz) * u.km / u.s

        return np.around(pos_mag, 3), np.around(v_mag, 3), m

    def all_particle_properties(self, particle_type=None, as_table=True):
        """
        Calculates distances from the origin and magnitude of the velocities
        for all particles (default) or a specied particle type.
        
        Kwargs: 
            particle_type (int):
                A subset of the data filtered by 1=DM, 2=disk, 3=bulge
            as_table (boolean): Return type. 
                If True, an astropy QTable with units. 
                If False, np.ndarrays for position and velocity

        Returns:
            QTable: 
            The full list, optionally with units, optionally filtered by type.
        """

        if particle_type is None:  # all types accepted
            dataset = self.data
        else:
            dataset = self.filter_by_type(particle_type)

        # mass:
        m = dataset['m'] * 1e10 * u.Msun

        # Pythagorean distance from origin:
        pos_xyz = np.array([dataset['x'], dataset['y'], dataset['z']])
        pos_mag = norm(pos_xyz, axis=0) * u.kpc

        # velocity:
        v_xyz = np.array([dataset['vx'], dataset['vy'], dataset['vz']])
        v_mag = norm(v_xyz, axis=0) * u.km / u.s

        if as_table:
            # construct and return a QTable with units
            t = QTable()
            t['type'] = dataset['type']
            t['m'] = np.around(m)
            t['pos'] = np.around(pos_mag, 3)
            t['v'] = np.around(v_mag, 3)
            return t
        else:
            return dataset['type'], dataset['m'], pos_mag, v_mag

    def component_count(self, particle_type=None):
        """
        Kwargs: particle_type (int):
                a subset of the data filtered by 1=DM, 2=disk, 3=bulge

        Returns: Quantity:
                The number of particles in the galaxy of this type
        """

        if particle_type is None:  # all types accepted
            dataset = self.data
        else:
            dataset = self.filter_by_type(particle_type)

        return len(dataset['m'])

    def all_component_counts(self):
        """
        Returns: list:
                The aggregate masses of particles of each type in the galaxy
                Ordered as [halo, disk, bulge]
        """

        return [self.component_count(ptype) for ptype in (1, 2, 3)]

    def component_mass(self, particle_type=None):
        """
        Kwargs: particle_type (int):
                a subset of the data filtered by 1=DM, 2=disk, 3=bulge

        Returns: Quantity:
                The aggregate mass of all particles in the galaxy of this type
        """

        if particle_type is None:  # all types accepted
            dataset = self.data
        else:
            dataset = self.filter_by_type(particle_type)

        return np.sum(dataset['m']) * 1e10 * u.Msun

    def all_component_masses(self):
        """
        Returns: list:
                The aggregate masses of particles of each type in the galaxy
        """

        return [self.component_mass(ptype) for ptype in (1, 2, 3)]

        # --------------------------------------------------------------------
    # Define some getters which may turn out to be useful, perhaps
    # --------------------------------------------------------------------

    def get_array(self):
        """
        Returns: all particle data in `np.ndarray` format

        Pretty superfluous in Python (which has no private class members)
        """

        return self.data

    def get_df(self):
        """
        Returns: data as pandas dataframe
        """

        return pd.DataFrame(self.data)

    def get_qtable(self):
        """
        Returns: data as astropy QTable, with units
        """

        t = QTable(self.data)

        # add appropriate units
        t['m'] = np.around(t['m'] * 1e10) * u.Msun
        for col in ('x', 'y', 'z'):
            t[col] *= u.kpc
        for col in ('vx', 'vy', 'vz'):
            t[col] *= (u.km / u.s)

        return t

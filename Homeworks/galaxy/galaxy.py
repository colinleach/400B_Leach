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
    A class to find, read and manipulate files for a single galaxy
    
    Needs to be initialized with a name, eg 'MW', 'M31'
    
    Snap number is optional, defaults to zero

    Path to the data file is optional, it will search several plausible defaults
    """

    def __init__(self, name, snap=0, datadir=None):
        "Initial setup"

        self.name = name
        self.snap = snap

        # We can probably make some assumptions about the filename:
        self.filename = f"{name}_{snap:03}.txt"
        self.filepath = self.get_filepath(datadir)
        self.read_file()

    def get_filepath(self, datadir):
        """
        Search for the required file and return a valid pathlib.Path object.
        
        This is pretty boring housekeeping code but may make things more resilient.
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
            return Path(datadir) # whoever called this clearly knew the layout

        # Now we aren't sure where the file is so may need to search for it
        # it may be different on different systems, e.g. nimoy vs laptop
        # try the search path [datadir, '.', '..', '../..', '~', '~/galaxydata']
        pwd = Path.cwd()
        home = Path.home()
        searchpath = [pwd, ]  
        searchpath += [pwd.parents[n] for n in range(4)]
        searchpath += [home, home / 'galaxydata']

        for p in searchpath:
            if has_file(p):
                return p # happy result, we got a valid path
            
        # Raise an error if not found
        pathstrings = f'{datadir}, ' if datadir is not None else ''
        pathstrings += ', '.join([sp.as_posix() for sp in searchpath])
        msg = f'{self.filename}: Not found in {pathstrings}'
        raise FileNotFoundError(msg)
        
    def read_file(self):
        """
        Read in a datafile in np.ndarray format, store in `self.data`
        
        Assumes path and filename are already set as instance parameters
        
        Returns: nothing
        """
        
        fullname = self.filepath / self.filename
        
        # get header data
        with open(fullname) as file:
            label, value = file.readline().split()
            self.time = float(value) * 10.0 * u.Myr
            label, value = file.readline().split()
            self.particle_count = int(value)
            
        # get the big array of values
        self.data = np.genfromtxt(fullname, dtype=None, names=True, skip_header=3)
        
    def filter_by_type(self, type, dataset=None):
        """
        Input:  particle type as integer, 1=DM, 2=disk, 3=bulge
                optionally, a starting dataset other than self.data
        
        Returns: subset data
        """
        
        if dataset is None:
            dataset = self.data
            
        index = np.where(dataset['type'] == type)
        return dataset[index]
        
    def single_particle_properties(self, particle_num=0, type=None):
        """
        Parameters:
            particle_num: zero-based index to an array of particles
            type: optionally, use a subset of the data filtered by
                    1=DM, 2=disk, 3=bulge
                    
            returns: 3-tuple of 
                        Euclidean distance from CoM (kpc), 
                        Euclidean velocity magnitude (km/s), 
                        particle mass (M_sun)
        """
        
        if type is None: # all types accepted
            particle = self.data[particle_num]
        else:
            particle = self.filter_by_type(type)[particle_num]
        
        # mass:
        m = particle['m'] * 1e10 * u.Msun
        
        # Euclidean distance from galactic CoM:
        pos_xyz = np.array([particle['x'], particle['y'], particle['y']])
        pos_mag = norm(pos_xyz) * u.kpc
        
        # velocity:
        v_xyz = np.array([particle['vx'], particle['vy'], particle['vy']])
        v_mag = norm(v_xyz) * u.km/ u.s
        
        return np.around(pos_mag, 3), np.around(v_mag, 3), m
    
    def all_particle_properties(self, type=None):
        """
        Similar to single_particle_properties, 
        except returns the full list 
        optionally filtered by type
        """
        
        if type is None: # all types accepted
            dataset = self.data
        else:
            dataset = self.filter_by_type(type)
            
        # mass:
        m = dataset['m'] * 1e10 * u.Msun
        
        # Pythagorean distance from galactic CoM:
        pos_xyz = np.array([dataset['x'], dataset['y'], dataset['y']])
        pos_mag = norm(pos_xyz, axis=0) * u.kpc
        
        # velocity:
        v_xyz = np.array([dataset['vx'], dataset['vy'], dataset['vy']])
        v_mag = norm(v_xyz, axis=0) * u.km/ u.s
        
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
        Returns np.ndarray format
        
        Pretty superfluous in Python (which has no private class members)
        """
        
        return self.data
    
    def get_df(self):
        "Returns pandas dataframe"
        
        return pd.DataFrame(self.data)
    
    def get_qtable(self):
        "returns astropy QTable with units"
        
        t = QTable(self.data)
        
        # add appropriate units
        t['m'] = np.around(t['m']*1e10) * u.Msun
        for col in ('x', 'y', 'z'):
            t[col] *= u.kpc
        for col in ('vx', 'vy', 'vz'):
            t[col] *= (u.km/u.s)
        
        return t
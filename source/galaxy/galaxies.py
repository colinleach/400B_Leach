# scientific package imports
import numpy as np
from numpy.linalg import norm
import pandas as pd

import astropy.units as u
from astropy.table import QTable

from galaxy.galaxy import Galaxy

class Galaxies():
    """
    Kwargs:
        names (iterable of str): 
            short names used in filename of type 'name_000.txt', eg 'MW', 'M31'.
        snaps (iterable of int): 
            Snap number, equivalent to time elapsed. Zero is starting conditions.
        datadir (str):
            Directory to search first for the required file. Optional, and a 
            default list of locations will be searched.

    Class attributes:
        path (`pathlib.Path` object): 
            directory (probably) containing the data files
        filenames (list of str):
            in `name_snap` format, something like 'MW_000' (no extension)
        galaxies (dict):
            key is filename, value is the corresponding Galaxy object

    A class to manipulate data for multiple galaxies.
    """

    def __init__(self, 
                names=('MW', 'M31', 'M33'),
                snaps=(0, 0, 0), 
                datadir=None):
        "Initial setup. Currently it calls read_file(), but this may change."

        self.names = names
        self.snaps = snaps
        self.path = datadir
        self.galaxies = {}
        self.filenames = []

    def read_data_files(self):
        """
        Attempts to create a Galaxy object for each name/snap combination
        set in `self.names` and `self.snaps`

        No return value.
        Sets self.galaxies, a dictionary keyed on name_snap
        """

        for name, snap in zip(self.names, self.snaps):
            # build the very important dictionary:
            key = f'{name}_{snap:03}' # e.g 'MW_000'
            self.galaxies[key] = Galaxy(name, snap, self.path)

            # bits of minor housekeeping:
            self.path = self.galaxies[key].path # may speed up next search
            self.filenames.append(key)

    def get_counts(self):
        """
        Args: none

        Returns: 
        """

        pass

    def get_masses(self):
        """
        Args: none

        Returns: 
        """

        pass

    def get_full_df(self):
        """
        Args: none

        Returns: 
            Concatenated pandas dataframe from all galaxies
            Includes 'name' and 'snap' columns
        """

        galaxies = []
        for i, gal_name in enumerate(self.filenames):
            g_df = self.galaxies[gal_name].all_particle_properties().to_pandas()
            g_df['name'] = self.names[i]
            g_df['snap'] = self.snaps[i]
            galaxies.append(g_df)
        return pd.concat(galaxies)
       
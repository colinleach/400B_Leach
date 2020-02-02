# scientific package imports
import numpy as np
from numpy.linalg import norm
import pandas as pd

import astropy.units as u
from astropy.table import QTable, Table

from galaxy import Galaxy

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
        "Initial setup."

        self.names = names
        self.snaps = snaps
        self.path = datadir
        self.galaxies = {}
        self.filenames = []
        self.read_data_files()

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
            self.path = self.galaxies[key].filepath # may speed up next search
            self.filenames.append(key)

    def get_counts_pivot(self):
        """
        Args: none

        Returns: 
        """

        gals_df = self.get_full_df()
        types = {1: '1 Halo', 2: '2 Disk', 3: '3 Bulge'}
        gals_df['typename'] = gals_df['type'].map(types)
        df_piv = pd.pivot_table(gals_df, values='m', index='name', 
            columns='typename', aggfunc='count', fill_value=0, margins=True)
            
        return df_piv

    def get_masses_pivot(self):
        """
        Args: none

        Returns: 
        """

        gals_df = self.get_full_df()
        types = {1: '1 Halo', 2: '2 Disk', 3: '3 Bulge'}
        gals_df['typename'] = gals_df['type'].map(types)
        df_piv = pd.pivot_table(gals_df, values='m', index='name', 
            columns='typename', aggfunc='sum', fill_value=0, margins=True)
            
        return df_piv

    def get_masses(self):
        """
        Args: none

        Returns: 
        """

        t = Table()
        t['Galaxy Name'] = self.names
        halo = []
        disk = []
        bulge = []
 
        for fname in self.filenames:
            h, d, b = self.galaxies[fname].all_component_masses()
            halo.append(h)
            disk.append(d)
            bulge.append(b)

        t['Halo Mass'] = halo
        t['Disk Mass'] = disk
        t['Bulge Mass'] = bulge
        
        lg_halo = sum(t['Halo Mass'])
        lg_disk = sum(t['Disk Mass'])
        lg_bulge = sum(t['Bulge Mass'])
        
        totals = ['Local Group', lg_halo, lg_disk, lg_bulge]
        t.add_row(totals)
        t['Total'] = t['Disk Mass'] + t['Bulge Mass'] + t['Halo Mass']
        t['f_bar'] = (t['Disk Mass'] + t['Bulge Mass']) / t['Total']

        return t


    def get_masses_df(self):
        """
        Args: none

        Returns: 
        """

        df = pd.DataFrame()
        # names = list(self.names) + ['Local Group',]
        df['Galaxy Name'] = list(self.names) + ['Local Group',]
        halo = []
        disk = []
        bulge = []
 
        for fname in self.filenames:
            h, d, b = self.galaxies[fname].all_component_masses()
            halo.append(h.value)
            disk.append(d.value)
            bulge.append(b.value)
        halo.append(sum(halo))
        disk.append(sum(disk))
        bulge.append(sum(bulge))
        

        df['Halo Mass (Msun*1e12)'] = np.around(np.array(halo) / 1e12, 3)
        df['Disk Mass (Msun*1e12)'] = np.around(np.array(disk) / 1e12, 3)
        df['Bulge Mass (Msun*1e12)'] = np.around(np.array(bulge) / 1e12, 3)
        
        lg_halo = sum(df['Halo Mass (Msun*1e12)'])
        lg_disk = sum(df['Disk Mass (Msun*1e12)'])
        lg_bulge = sum(df['Bulge Mass (Msun*1e12)'])
        
        totals = ['Local Group', lg_halo, lg_disk, lg_bulge]
        # df.add_row(totals)
        df['Total'] = df['Disk Mass (Msun*1e12)'] \
            + df['Bulge Mass (Msun*1e12)']  + df['Halo Mass (Msun*1e12)']
        df['f_bar'] = np.around((df['Disk Mass (Msun*1e12)'] + \
            df['Bulge Mass (Msun*1e12)']) / df['Total'], 3)

        return df


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
       
# scientific package imports
import numpy as np
import pandas as pd
from numpy.linalg import norm

import astropy.units as u
from astropy.table import QTable #, Table

from galaxy.galaxy import Galaxy
from galaxy.centerofmass import CenterOfMass


class Galaxies():
    """
    A class to manipulate data for multiple galaxies.

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
        Sets `self.galaxies`, a dictionary keyed on `name_snap`
        """

        for name, snap in zip(self.names, self.snaps):
            # build the very important dictionary:
            key = f'{name}_{snap:03}'  # e.g 'MW_000'
            self.galaxies[key] = Galaxy(name, snap, self.path)

            # bits of minor housekeeping:
            self.path = self.galaxies[key].filepath  # may speed up next search
            self.filenames.append(key)

    def get_pivot(self, aggfunc, values='m'):
        """
        Generic method to make a pandas pivot table from the 9 combinations of 
        galaxy and particle type.

        Args:
            aggfunc (str): 'count', 'sum', etc as aggregation method
            values (str): column name to aggregate

        Returns: pandas dataframe
        """

        # stack all galaxies in a pandas dataframe
        gals_df = self.get_full_df()

        # create some better column names
        types = {1: '1 Halo', 2: '2 Disk', 3: '3 Bulge'}
        gals_df['typename'] = gals_df['type'].map(types)

        # get pandas to do most of the work
        df_piv = pd.pivot_table(gals_df, values=values,
                                index='name', columns='typename',
                                aggfunc=aggfunc, fill_value=0, margins=True)

        return df_piv

    def get_counts_pivot(self):
        """
        Pivots on `count('m)`.

        Returns: pandas dataframe
        """

        return self.get_pivot('count')

    def get_masses_pivot(self):
        """
        Pivots on `sum('m)`.
        
        Returns: pandas dataframe
        """

        return self.get_pivot('sum')

    def get_full_df(self):
        """
        Combined data for all input files.

        Returns:
            Concatenated pandas dataframe from all galaxies
            Includes 'name' and 'snap' columns
        """

        galaxies = []
        for i, gal_name in enumerate(self.filenames):
            g_df = self.galaxies[gal_name].all_particle_properties(
                ).to_pandas()
            g_df['name'] = self.names[i]
            g_df['snap'] = self.snaps[i]
            galaxies.append(g_df)
        return pd.concat(galaxies)

    def get_coms(self, tolerance=0.1, ptypes=(1,2,3)):
        """
        Center of Mass determination for all galaxies. 
        Defaults to all particle types, but `ptypes=(2,)` may be more useful.

        Args:
            tolerance (float): convergence criterion (kpc)

        Returns:
            QTable with COM positions and velocities
            colnames: ['name', 'ptype', 'x', 'y', 'z', 'vx', 'vy', 'vz', 'R', 'V']
        """

        # gather the position/velocity data for all galaxies and particle types
        vals = []

        for name in self.filenames:
            for ptype in ptypes:
                g = self.galaxies[name]
                try:
                    com = CenterOfMass(g, ptype)
                except ValueError as err:
                    xyz_com = (np.NaN, np.NaN, np.NaN) * u.kpc
                    vxyz_com = (np.NaN, np.NaN, np.NaN) * u.km / u.s
                else:
                    xyz_com= com.com_p(tolerance)
                    vxyz_com= com.com_v(xyz_com)
                row = [name, g.type2name(ptype), xyz_com, vxyz_com]
                vals.append(row)

        # The QTable constructor is fussy and tends to surprise
        # first transpose the list of lists
        # ref: https://stackoverflow.com/questions/6473679/transpose-list-of-lists
        vals_transposed = list(map(list, zip(*vals)))

        # We can build a QTable with 3-vector entries
        names = ('name', 'ptype', 'xyz', 'vxyz')
        qt = QTable(vals_transposed, names=names)

        # Construct a better QTable, one value per column
        # I feel there ought to be an easier way than this?
        qt2 = QTable()
        qt2['name'] = qt['name']
        qt2['ptype'] = qt['ptype']
        qt2['x'] = [x for x, y, z in qt['xyz']]
        qt2['y'] = [y for x, y, z in qt['xyz']]
        qt2['z'] = [z for x, y, z in qt['xyz']]
        qt2['vx'] = [x for x, y, z in qt['vxyz']]
        qt2['vy'] = [y for x, y, z in qt['vxyz']]
        qt2['vz'] = [z for x, y, z in qt['vxyz']]

        # Add columns for distance from origin and velocity magnitude
        # Using `norm` failed with a strange dtype conversion error
        # so do it the old-school way
        qt2['R'] = np.round(np.sqrt(qt2['x']**2 + qt2['y']**2 + qt2['z']**2), 2)
        qt2['V'] = np.round(np.sqrt(qt2['vx']**2 + qt2['vy']**2 + qt2['vz']**2), 2)

        return qt2

    def separations(self, g1, g2):
        """
        Position and velocity of galaxy g2 COM relative to g1 COM. 
        Uses only disk particles for the COM determination.

        Args:
            g1, g2 (str): galaxies matching entries in self.filenames

        Returns:
            Dictionary containing relative position, distance, velocities in
            Cartesian and radial coordinates
        """

        results = {}
        com1 = CenterOfMass(self.galaxies[g1], 2)
        com2 = CenterOfMass(self.galaxies[g2], 2)
        com1_p = com1.com_p()
        com2_p = com2.com_p()
        com1_v = com1.com_v(com1_p)
        com2_v = com2.com_v(com2_p)
        
        results['pos_xyz'] = com2_p - com1_p
        results['vel_xyz'] = com2_v - com1_v
        results['r'] = np.round(norm(results['pos_xyz']), 2)
        results['r_hat'] = np.round(results['pos_xyz'] / results['r'], 2) # unit vector
        results['vel_mag'] = np.round(norm(results['vel_xyz']), 2)
        results['v_radial'] = np.round(np.dot(results['r_hat'], results['vel_xyz']), 2)
        results['v_tangential'] = np.round(np.cross(results['r_hat'], results['vel_xyz']), 2)
        results['v_tan_mag'] = np.round(norm(results['v_tangential']), 2)

        return results

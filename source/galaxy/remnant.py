# standard Python imports
from pathlib import Path

# scientific package imports
import numpy as np
from numpy.linalg import norm
import astropy.units as u

from galaxy.galaxy import Galaxy


class Remnant(Galaxy):
    """
    A class to work with the MW-M31 post-merger remant.

    Args:
        snap (int):
            Snap number, equivalent to time elapsed. 
            Defaults to the last timepoint.
        datadir (str):
            Directory to search first for the required file. Optional, and a
            default list of locations will be searched.
        usesql (bool):
            If True, data will be taken from a PostgreSQL database instead of
            text files.
        stride (int):
            Optional. For stride=n, get every nth row in the table.
            Only valid with usesql=True.

    Class attributes:
        data (np.ndarray):
            type, mass, position_xyz, velocity_xyz for each particle
    """

    def __init__(self, snap=801, datadir=None, usesql=False, stride=1):
        "Initial setup. Currently it calls read_file(), but this may change."

        self.snap = snap

        if usesql:
            self.read_db(stride)
        else:
            raise NotImplementedError

    def read_db(self, stride):
        """
        Get relevant data from a PostgreSQL database and format it to be 
        identical to that read from test files.

        Args:
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
        sql_t = f"SELECT time FROM simdata WHERE galname in ('MW', 'M31')"
        sql_t += f" and snap={self.snap} LIMIT 1"
        cur.execute(sql_t)
        time = cur.fetchone()
        try:
            self.time = time[0] * u.Myr
        except TypeError:
            print(self.name, self.snap, ptype)

        # set the bulk of the data
        colheads = ','.join(['galname','type','m','x','y','z','vx','vy','vz'])
        if stride > 1:
            sql_d = f"SELECT {colheads}, ROW_NUMBER() OVER () as rn from simdata"
        else:
            sql_d = f"SELECT {colheads} from simdata"
        sql_d += f"  where galname in ('MW', 'M31') and snap={self.snap}"
        sql_d += f" and type in (2,3)"
        sql_d += " ORDER BY galname, pnum"
        if stride > 1:
            sql_d = f"SELECT {colheads} from ( {sql_d} ) as t where rn % {stride} = 0" 

        dtype=[('galname', 'U3'), ('type', 'uint8'), ('m', '<f4'), 
                ('x', '<f4'), ('y', '<f4'), ('z', '<f4'), 
                ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        cur.execute(sql_d)
        self.data = np.array(cur.fetchall(), dtype=dtype)
        self.particle_count = len(self.data)
        
    def xyz(self):
        """
        Convenience method to get positions as a np.array of shape (3,N)
        """

        return np.array([self.data[xi] for xi in ('x','y','z')])

    def vxyz(self):
        """
        Convenience method to get velocities as a np.array of shape (3,N)
        """
        
        return np.array([self.data[vxi] for vxi in ('vx','vy','vz')])

        
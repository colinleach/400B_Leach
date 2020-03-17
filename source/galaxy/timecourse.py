# import modules
from pathlib import Path
import numpy as np
from numpy.linalg import norm

import astropy.units as u

from galaxy.galaxy import Galaxy
from galaxy.galaxies import Galaxies
from galaxy.centerofmass import CenterOfMass
from galaxy.db import DB

class TimeCourse():

    def __init__(self, datadir='.', usesql=False):
        self.datadir = datadir
        self.usesql = usesql

    def write_com_ang_mom(self, galname, start=0, end=801, n=5, show_progress=True):
        """
        Function that loops over all the desired snapshots to compute the COM pos and vel as a 
        function of time.

        inputs:
            galname (str):
                'MW', 'M31' or 'M33'
            start, end (int):
                first and last snap numbers to include
            n (int):
                stride length for the sequence
            datadir (str):
                path to the input data
            show_progress (bool):
                prints each snap number as it is processed
            
        returns: 
            Two text files saved to disk. 
        """
        
        # compose the filenames for output
        com_fileout = f'./com_{galname}.txt'
        angmom_fileout = f'./angmom_{galname}.txt'
        
        #  set tolerance and vol_dec for calculating COM_P in CenterOfMass
        # for M33 that is stripped more, use different values for vol_dec
        if galname == 'M33':
            delta = 0.1
            vol_dec = 4
        else:
            delta = 0.2
            vol_dec = 2
        
         # generate the snapshot id sequence 
        snap_ids = np.arange(start, end+1, n)
        
        # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
        orbit = np.zeros((len(snap_ids), 7))
        
        # initialize the array for angular momentum info:
        #  t, (x, y, z)_hat, magnitude of L
        angmom = np.zeros((len(snap_ids), 5))
        
        if show_progress:
            print(galname)
        
        for  i, snap in enumerate(snap_ids): # loop over files
            gal = Galaxy(galname, snap, datadir=self.datadir, usesql=self.usesql)
            
            # Initialize an instance of CenterOfMass class, using disk particles
            com = CenterOfMass(gal)

            # Store the COM pos and vel. Remember that now COM_P required VolDec
            com_p = com.com_p(delta, vol_dec)
            com_v = com.com_v(com_p)

            # store the angular momentum as unit vector and magnitude
            L, _, _ = com.angular_momentum()
            L_mag = norm(L)
            L_hat = L / L_mag
        
            # store the time, pos, vel in ith element of the orbit array, without units
            orbit[i] = gal.time.value/1000, *tuple(com_p.value), *tuple(com_v.value)
            
           # store the time, L_hat, L_mag in ith element of the angular momentum array
            angmom[i] = gal.time.value/1000, *tuple(L_hat), L_mag
            
            # print snap_id to see the progress
            if show_progress:
                print(snap, end=' ')
        print('\nDone')
            
        # write the data to 2 files
        # we do this because we don't want to have to repeat this process 
        # this code should only have to be called once per galaxy.
        np.savetxt(com_fileout, orbit, fmt = "%11.3f"*7, comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

        np.savetxt(angmom_fileout, angmom, fmt = "%11.3f"*4 + "%11.3e", comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'x_hat', 'y_hat', 'z_hat', 'L_mag'))

    def write_total_com(self, start=0, end=801, n=5, show_progress=True):
        """
        """

        fileout = './total_com.txt'
        snap_ids = np.arange(start, end+1, n)
        coms = np.zeros((len(snap_ids), 7))

        for  i, snap in enumerate(snap_ids): # loop over files
            gals = Galaxies(snaps=[snap]*3, datadir=self.datadir, usesql=self.usesql)
            full_com_p, full_com_v = gals.total_com()
            t = gals.time.value/1000
            coms[i] = t, *tuple(full_com_p), *tuple(full_com_v)

            if show_progress:
                print(snap, end=' ')
        print('\nDone')

        np.savetxt(fileout, coms, fmt = "%11.3f"*7, comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))


    def read_file(self, fullname):
        """
        General method for file input. Note that the format is for summary files,
        (one line per snap), not the raw per-particle files.
        """

        with open(fullname) as file:
            data = np.genfromtxt(fullname, dtype=None, names=True, skip_header=0)
        return data
 

    def read_com_file(self, galaxy, datadir='.'):
        """
        Get CoM summary from file.

        Args:
            galaxy (str): 
                'MW', 'M31', 'M33'
            datadir (str):
                path to file

        Returns:
            np.array with 802 rows, one per snap
        """

        filename = f'com_{galaxy}.txt'
        fullname = Path(datadir) / filename

        return self.read_file(fullname)
    
    def read_angmom_file(self, galaxy, datadir='.'):
        """
        Get CoM summary from file.

        Args:
            galaxy (str): 
                'MW', 'M31', 'M33'
            datadir (str):
                path to file

        Returns:
            np.array with 802 rows, one per snap
        """

        filename = f'angmom_{galaxy}.txt'
        fullname = Path(datadir) / filename

        return self.read_file(fullname)
    
    def read_total_com_file(self, galaxy, datadir='.'):
        """
        Get CoM summary from file.

        Args:
            galaxy (str): 
                'MW', 'M31', 'M33'
            datadir (str):
                path to file

        Returns:
            np.array with 802 rows, one per snap
        """

        filename = f'com_{galaxy}.txt'
        fullname = Path(datadir) / filename

        return self.read_file(fullname)
    
    def write_db_tables(self, datadir='.'):
        """
        Adds data to the `centerofmass`, `angmom` and `totalcom` tables in the 
        `galaxy` database
        """

        filepath = Path(datadir)

        db = DB()
        cur = db.get_cursor()

        # CoM data
        colheads = ','.join(['gal','snap','t','x','y','z','vx','vy','vz'])
        query = f"""
            INSERT INTO centerofmass( {colheads} ) 
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s,%s)
            ON CONFLICT DO NOTHING
            """

        for gname in ('MW','M31','M33'):
            filename = f'com_{gname}.txt'
            fullname = filepath / filename
            data = self.read_file(fullname)
               
            for snap, d in enumerate(data):
                rec = [gname, snap,] + list(d)
                cur.execute(query, rec)

        # angular momentum data
        colheads = ','.join(['gal','snap','t','x_hat','y_hat','z_hat','l_mag'])
        query = f"""
            INSERT INTO angmom( {colheads} ) 
            VALUES (%s,%s,%s,%s,%s,%s,%s)
            ON CONFLICT DO NOTHING
            """

        for gname in ('MW','M31','M33'):
            filename = f'angmom_{gname}.txt'
            fullname = filepath / filename
            data = self.read_file(fullname)
                
            for snap, d in enumerate(data):
                rec = [gname, snap,] + list(d)
                cur.execute(query, rec)

        # total CoM data
        colheads = ','.join(['snap','t','x','y','z','vx','vy','vz'])
        query = f"""
            INSERT INTO totalcom( {colheads} ) 
            VALUES (%s,%s,%s,%s,%s,%s,%s,%s)
            ON CONFLICT DO NOTHING
            """

        filename = 'total_com.txt'
        fullname = filepath / filename
        data = self.read_file(fullname)
            
        for snap, d in enumerate(data):
            rec = [snap,] + list(d)
            cur.execute(query, rec)

    def read_com_db(self, galaxy=None, snaprange=(0,801)):
        """
        """

        colheads = ','.join(['gal','snap','t','x','y','z','vx','vy','vz'])
        query = f"""
                SELECT {colheads} FROM centerofmass 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                ORDER BY snap
                """
        if galaxy is not None:
            query += f" AND gal={galaxy}"

        db = DB()
        result = db.run_query(query)
        dtype=[('gal', 'U3'), ('snap', 'u2'), ('t', '<f4'), ('x', '<f4'), ('y', '<f4'), 
                ('z', '<f4'), ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        return np.array(result, dtype=dtype)
        
    def read_angmom_db(self, galaxy=None, snaprange=(0,801)):
        """
        """

        colheads = ','.join(['gal','snap','t','x_hat','y_hat','z_hat','l_mag'])
        query = f"""
                SELECT {colheads} FROM angmom 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                ORDER BY snap
                """
        if galaxy is not None:
            query += f" AND gal={galaxy}"

        db = DB()
        result = db.run_query(query)
        dtype=[('gal', 'U3'), ('snap', 'u2'), ('t', '<f4'), ('x_hat', '<f4'), ('y_hat', '<f4'), 
                ('z_hat', '<f4'), ('l_mag', '<f4')]

        return np.array(result, dtype=dtype)

    def read_total_com_db(self, snaprange=(0,801)):
        """
        """

        colheads = ','.join(['snap','t','x','y','z','vx','vy','vz'])
        query = f"""
                SELECT {colheads} FROM totalcom 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                ORDER BY snap
                """

        db = DB()
        result = db.run_query(query)
        dtype=[('snap', 'u2'), ('t', '<f4'), ('x', '<f4'), ('y', '<f4'), 
                ('z', '<f4'), ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        return np.array(result, dtype=dtype)
        
       
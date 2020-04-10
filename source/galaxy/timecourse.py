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
    """
    A commection of methods for generating, reading and writing summary data for 
    parameters that change over the timecourse of the simulation.

    These fall into a few groups:
        `write_xxx()` : 
            Methods that loop through the raw data, calculate parameters and write 
            the results to file. Can be slow (hours) to run but Only run once. 
            See the `data` folder for the resulting files, one line per snap.

        `read_xxx_file(`)` :
            Read in the summary files and return a numpy array. These rely on the
            generic `read_file()` method.

        `read_xxx_db()` :
            Get the summary data from postgres instead of a text file. 
            The returned array should be identical to the `read_xxx_file()` group.

        `write_db_tables()` :
            Read a text file, insert the contents to a postgres table.

        `get_one_com()` :
            Convenience method to return a single CoM position.
    """

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

    def write_total_com(self, start=0, end=801, n=1, show_progress=True):
        """
        Function that loops over all the desired snapshots to compute the overall COM 
        pos and vel as a function of time. Uses all particles in all galaxies.

        inputs:
            start, end (int):
                first and last snap numbers to include
            n (int):
                stride length for the sequence
            show_progress (bool):
                prints each snap number as it is processed
            
        output: 
            Text file saved to disk. 
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

    def write_total_angmom(self, start=0, end=801, n=1, show_progress=True):
        """
        Function that loops over all the desired snapshots to compute the overall 
        angular momentum as a function of time. Uses all particles in all galaxies.

        inputs:
            start, end (int):
                first and last snap numbers to include
            n (int):
                stride length for the sequence
            show_progress (bool):
                prints each snap number as it is processed
            
        output: 
            Text file saved to disk. 
        """

        fileout = './total_angmom.txt'
        snap_ids = np.arange(start, end+1, n)
        angmoms = np.zeros((len(snap_ids), 4))

        com = self.read_total_com_db([start, end])
        # print(com.shape, com['x'].shape)
        com_xyz = np.array([com[xi] for xi in ('x','y','z')]).T
        # print(com_xyz, com_xyz.shape)

        for  i, snap in enumerate(snap_ids): # loop over files
            gals = Galaxies(snaps=[snap]*3, datadir=self.datadir, usesql=self.usesql)
            L = gals.total_angmom(com_xyz[i])
            t = gals.time.value/1000
            angmoms[i] = t, *tuple(L)

            if show_progress:
                print(snap, end=' ')
        print('\nDone')

        np.savetxt(fileout, angmoms, fmt = "%11.3f"*4, comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'Lx', 'Ly', 'Lz'))

    def write_vel_disp(self, galname, start=0, end=801, n=1, show_progress=True):
        """
        Function that loops over all the desired snapshots to compute the veocity dispersion
        sigma as a function of time.

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
            Text file saved to disk. 
        """

        # compose the filename for output
        # fileout = f'./sigma_{galname}.txt'
        fileout = f'./sigma_{galname}_2.txt'

         # generate the snapshot id sequence 
        snap_ids = np.arange(start, end+1, n)
        
        # initialize the array for orbital info: t, sigma
        sigmas = np.zeros((len(snap_ids), 2))
        
        if show_progress:
            print(galname)
        
        for  i, snap in enumerate(snap_ids): # loop over files
            gal = Galaxy(galname, snap, datadir=self.datadir, usesql=self.usesql)
            
            # Initialize an instance of CenterOfMass class, using disk particles
            com = CenterOfMass(gal)

            # COM pos and vel
            com_xyz, com_vxyz = self.get_one_com(galname, snap)  

            gal_xyzD, gal_vxyzD = com.center_com(com_xyz, com_vxyz)

            # determine the rotated velocity vectors
            _, vn = com.rotate_frame(com_p=com_xyz, com_v=com_vxyz)
            
            # calculate velocity dispersion
            v_radial = vn[0] # just x-component
            v_mean = np.mean(v_radial)
            sigmas[i] = gal.time.value/1000, np.std(v_radial - v_mean)

            # print snap_id to see the progress
            if show_progress:
                print(snap, end=' ')
        print('\nDone')
            
        # write the data to file
        # we do this because we don't want to have to repeat this process 
        # this code should only have to be called once per galaxy.
        np.savetxt(fileout, sigmas, fmt = "%11.3f"*2, comments='#',
                header="{:>10s}{:>11s}"\
                        .format('t', 'sigma'))

    def write_LG_normal(self, start=0, end=801):
        """
        Calculates the normal to a plane containing the three galaxy CoMs.

        Args:
            start, end (int):
                first and last snap numbers to include
            
        output: 
            Text file saved to disk. 
        """

        # get all the CoM info from postgres
        MW_data = self.read_com_db('MW')
        M31_data = self.read_com_db('M31')
        M33_data = self.read_com_db('M33')
        
        # pull out just the 3 columns giving position
        MW_coms = np.array([MW_data[xi] for xi in ('x','y','z')])
        M31_coms = np.array([M31_data[xi] for xi in ('x','y','z')])
        M33_coms = np.array([M33_data[xi] for xi in ('x','y','z')])

        # define 2 vectors that lie in the plane
        M31_MW = MW_coms - M31_coms
        M31_M33 = M33_coms - M31_coms

        # the normal we want comes from the vector cross product
        normals = np.cross(M31_MW, M31_M33, axis=0)
        normals /= norm(normals, axis=0)

        output = np.concatenate((MW_data['t'][:,np.newaxis], normals.T), axis=1)
        print(output.shape)
            
        # compose the filename for output
        fileout = './normals.txt'

        # write the data to file
        # we do this because we don't want to have to repeat this process 
        # this code should only have to be called once per galaxy.
        np.savetxt(fileout, output, fmt = "%11.3f"*4, comments='#',
                header="{:>10s}{:>11s}{:>11s}{:>11s}"\
                        .format('t', 'x_hat', 'y_hat', 'z_hat'))
              

    def read_file(self, fullname):
        """
        General method for file input. Note that the format is for summary files,
        (one line per snap), not the raw per-particle files.
        """

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
    
    def read_total_com_file(self, datadir='.'):
        """
        Get CoM summary from file.

        Args:
            datadir (str):
                path to file

        Returns:
            np.array with 802 rows, one per snap
        """

        filename = 'total_com.txt'
        fullname = Path(datadir) / filename

        return self.read_file(fullname)
    
    def read_normals_file(self, datadir='.'):
        """
        Get normals to plane containing 3 galaxy CoMs from file.

        Args:
            datadir (str):
                path to file

        Returns:
            np.array with 802 rows, one per snap
        """

        filename = 'normals.txt'
        fullname = Path(datadir) / filename

        return self.read_file(fullname)

    def write_db_tables(self, datadir='.', do_com=False, do_angmom=False, 
                        do_totalcom=False, do_totalangmom=False, do_normals=False,
                        do_sigmas=False):
        """
        Adds data to the various tables in the `galaxy` database
        """

        filepath = Path(datadir)

        db = DB()
        cur = db.get_cursor()

        # CoM data
        if do_com:
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
        if do_angmom:
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
        if do_totalcom:
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

        # total angular momentum data
        if do_totalangmom:
            colheads = ','.join(['snap','t','Lx','Ly','Lz'])
            query = f"""
                INSERT INTO totalangmom( {colheads} ) 
                VALUES (%s,%s,%s,%s,%s)
                ON CONFLICT DO NOTHING
                """

            filename = 'total_angmom.txt'
            fullname = filepath / filename
            data = self.read_file(fullname)
                
            for snap, d in enumerate(data):
                rec = [snap,] + list(d)
                cur.execute(query, rec)

        # 3-galaxy normals data
        if do_normals:
            colheads = ','.join(['snap','t','x_hat','y_hat','z_hat'])
            query = f"""
                INSERT INTO normals( {colheads} ) 
                VALUES (%s,%s,%s,%s,%s)
                ON CONFLICT DO NOTHING
                """

            filename = 'normals.txt'
            fullname = filepath / filename
            data = self.read_file(fullname)
                
            for snap, d in enumerate(data):
                rec = [snap,] + list(d)
                cur.execute(query, rec)

        # velocity dispersions
        if do_sigmas:
            colheads = ','.join(['gal','snap','t','sigma'])
            query = f"""
                INSERT INTO sigmas( {colheads} ) 
                VALUES (%s,%s,%s,%s)
                ON CONFLICT DO NOTHING
                """

            for gname in ('MW','M31','M33'):
                filename = f'sigma_{gname}.txt'
                fullname = filepath / filename
                data = self.read_file(fullname)
                    
                for snap, d in enumerate(data):
                    rec = [gname, snap,] + list(d)
                    cur.execute(query, rec)

    def read_com_db(self, galaxy=None, snaprange=(0,801)):
        """
        Retrieves CoM positions from postgres for a range of snaps.

        Args:
            galaxy (str):
                Optional, defaults to all. Can be 'MW', 'M31' , 'M33'
            snaprange (pair of ints):
                Optional, defaults to all. First and last snap to include.
                This is NOT the [first, last+1] convention of Python.
       """

        colheads = ','.join(['gal','snap','t','x','y','z','vx','vy','vz'])
        query = f"""
                SELECT {colheads} FROM centerofmass 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                """
        if galaxy is None:
            query += " ORDER BY snap"
        else:
            query += f" AND gal='{galaxy}' ORDER BY galaxy, snap"

        db = DB()
        result = db.run_query(query)
        dtype=[('gal', 'U3'), ('snap', 'u2'), ('t', '<f4'), ('x', '<f4'), ('y', '<f4'), 
                ('z', '<f4'), ('vx', '<f4'), ('vy', '<f4'), ('vz', '<f4')]

        return np.array(result, dtype=dtype)
        
    def read_angmom_db(self, galaxy=None, snaprange=(0,801)):
        """
        Retrieves disk angular momentum from postgres for a range of snaps.

        Args:
            galaxy (str):
                Optional, defaults to all. Can be 'MW, 'M31 , 'M33'
            snaprange (pair of ints):
                Optional, defaults to all. First and last snap to include.
                This is NOT the [first, last+1] convention of Python.
       """

        colheads = ','.join(['gal','snap','t','x_hat','y_hat','z_hat','l_mag'])
        query = f"""
                SELECT {colheads} FROM angmom 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
               """
        if galaxy is not None:
            query += f" AND gal='{galaxy}'"
        query += " ORDER BY snap"

        db = DB()
        result = db.run_query(query)
        dtype=[('gal', 'U3'), ('snap', 'u2'), ('t', '<f4'), ('x_hat', '<f4'), ('y_hat', '<f4'), 
                ('z_hat', '<f4'), ('l_mag', '<f4')]

        return np.array(result, dtype=dtype)

    def read_total_com_db(self, snaprange=(0,801)):
        """
        Retrieves total CoM positions from postgres for a range of snaps.

        Args:
            snaprange (pair of ints):
                Optional, defaults to all. First and last snap to include.
                This is NOT the [first, last+1] convention of Python.
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
        
    def get_one_com(self, gal, snap)       :
        """
        Gets a CoM from postgres for the specified galaxy and snap.

        Args:
            gal (str): 
                Can be 'MW, 'M31 , 'M33'
            snap (int):
                The timepoint.
        """

        com = self.read_com_db(gal, (snap,snap))
        xyz = np.array([com[xi] for xi in ('x','y','z')])
        vxyz = np.array([com[vxi] for vxi in ('vx','vy','vz')])

        # need to massage the shapes a bit to match what c.com_p() returns
        return xyz.T[0], vxyz.T[0]
        
    def read_total_angmom_db(self, snaprange=(0,801)):
        """
        Gets the total angular momentum of the 3-galaxy system. In practice, 
        this turns out to be near-zero at all timepoints and can be ignored
        in future work.
        """

        colheads = ','.join(['snap','t','Lx','Ly','Lz'])
        query = f"""
                SELECT {colheads} FROM totalangmom 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                ORDER BY snap
               """
 
        db = DB()
        result = db.run_query(query)
        dtype=[('snap', 'u2'), ('t', '<f4'), ('Lx', '<f4'), ('Ly', '<f4'), 
                ('Lz', '<f4')]

        return np.array(result, dtype=dtype)

    def read_normals_db(self, snaprange=(0,801)):
        """
        Gets the normals to the 3-galaxy plane.
        """

        colheads = ','.join(['snap','t','x_hat','y_hat','z_hat'])
        query = f"""
                SELECT {colheads} FROM normals 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                ORDER BY snap
               """
 
        db = DB()
        result = db.run_query(query)
        dtype=[('snap', 'u2'), ('t', '<f4'), ('x_hat', '<f4'), ('y_hat', '<f4'), 
                ('z_hat', '<f4')]

        return np.array(result, dtype=dtype)

    def read_sigmas_db(self, galaxy=None, snaprange=(0,801)):
        """
        Gets the velocity dispersions (km/s) for one galaxy at a range of snaps.
        """

        colheads = ','.join(['gal','snap','t','sigma'])
        query = f"""
                SELECT {colheads} FROM sigmas 
                WHERE snap BETWEEN {snaprange[0]} AND {snaprange[1]}
                """
        if galaxy is not None:
            query += f" AND gal='{galaxy}'"
        query += " ORDER BY snap"

        db = DB()
        result = db.run_query(query)
        dtype=[('gal', 'U3'), ('snap', 'u2'), ('t', '<f4'), ('sigma', '<f4')]

        return np.array(result, dtype=dtype)


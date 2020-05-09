# import modules
import numpy as np
from numpy.linalg import norm

import astropy.units as u
import astropy.table as tbl

from galaxy.galaxy import Galaxy

class CenterOfMass:
    """
    Class to define COM position and velocity properties of a given galaxy 
    and simulation snapshot

    Args:
        gal (Galaxy object):
            The desired galaxy/snap to operate on
        ptype (int):
            for particles, 1=DM/halo, 2=disk, 3=bulge

    Throws:
        ValueError, if there are no particles of this type in this galaxy
        (typically, halo particles in M33)
    """
    
    
    def __init__(self, gal, ptype=2):
    # Initialize the instance of this Class with the following properties:
    
        # subset with our particle type
        if ptype is None:
            self.data = gal.data
        else:
            self.data = gal.filter_by_type(ptype)     
            if self.data.shape[0] == 0:
                raise ValueError(f'No particles of type {ptype} in {gal.filename}')

        #create an array to store indexes of particles of desired Ptype                                
        # self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type
        self.m = self.data['m']
        self.xyz = np.array([self.data[col] for col in ('x', 'y', 'z')])
        self.vxyz = np.array([self.data[col] for col in ('vx', 'vy', 'vz')])


    def com_define(self, xyz, m):
        """
        Function to compute the center of mass position or velocity generically

        Args: 
            xyz (array with shape (3, N)):
                (x, y, z) positions or velocities
            m (1-D array):
                particle masses 

        Returns: 
            3-element array, the center of mass coordinates
        """

        return np.sum(xyz * m, axis=1) / np.sum(m)
    
    
    def com_p(self, delta=0.1, vol_dec=2.0):
        """
        Function to specifically return the center of mass position and velocity    .

        Kwargs:                                                                                                           
            delta (tolerance)                                                                                         
        Returns: 
            One 3-vector, coordinates of the center of mass position (kpc)   
        """                                                          

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling com_define                                                   
        xyz_com = self.com_define(self.xyz, self.m)
        # compute the magnitude of the COM position vector.
        RCOM = norm(xyz_com)

       # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        xyz_new = self.xyz - xyz_com[:, np.newaxis]
        r_new = norm(xyz_new, axis=0)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at reduced radius                                                         
        r_max = max(r_new)/vol_dec
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        change = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (change > delta):
            # select all particles within the reduced radius 
            # (starting from original x,y,z, m)
            index2 = np.nonzero(r_new < r_max)[0]
            xyz2 = self.xyz[:, index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            xyz_com2 = self.com_define(xyz2, m2)
            # compute the new 3D COM position
            r_com2 = norm(xyz_com2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            change = np.abs(RCOM - r_com2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            r_max = r_max/vol_dec
            # check this.                                                                                              
            # print ("RMAX", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            xyz_new = self.xyz - xyz_com2[:, np.newaxis]
            r_new = norm(xyz_new, axis=0)

            # set the center of mass positions to the refined values                                                   
            xyz_com = xyz_com2
            RCOM = r_com2

        # set the correct units using astropy and round all values
        # and then return the COM position vector
        return np.round(xyz_com, 2) * u.kpc
        
    def com_v(self, xyz_com):
        """
        Center of Mass velocity

        Args: X, Y, Z positions of the COM (no units)

        Returns: 3-Vector of COM velocities
        """

        try:
            xyz_com = xyz_com.value # in case units passed in
        except AttributeError:
            pass 
        
        # the max distance from the center that we will use to determine the CoM velocity                   
        rv_max = 15.0 # implicit u.kpc

        # determine the position of all particles relative to the center of mass position
        xyz_v = self.xyz - xyz_com[:, np.newaxis]
        r_v = norm(xyz_v, axis=0)
        
        # determine the index for those particles within the max radius
        index_v = np.where(r_v < rv_max)[0]

        # determine the velocity and mass of those particles within the max radius
        vxyz_new = self.vxyz[:, index_v]
        m_new = self.m[index_v]  
        
        # compute the center of mass velocity using those particles
        vxyz_com = self.com_define(vxyz_new, m_new)

        # return the COM vector                                                                                        
        return np.round(vxyz_com, 2) * u.km / u.s

    def center_com(self, com_p=None, com_v=None):
        """
        Positions and velocities of disk particles relative to the CoM

        Returns : two (3, N) arrays
            CoM-centric position and velocity
        """

        if com_p is None:
            com_p = self.com_p(0.1).value
        if com_v is None:
            com_v = self.com_v(com_p).value

        # Determine positions of disk particles relative to COM 
        xyzD = self.xyz - com_p[:, np.newaxis]

        # Determine velocities of disk particles relative to COM motion
        vxyzD = self.vxyz - com_v[:, np.newaxis]

        return xyzD, vxyzD


    def angular_momentum(self, com_p=None, com_v=None, r_lim=None):
        """
        Returns: 
            L : 3-vector as array
                The (x,y,x) components of the angular momentum vector about the CoM,
                summed over all disk particles
            pos, v : arrays with shape (3, N)
                Position and velocity for each particle
            r_lim : float
                Radius to include in calculation (implicit kpc, no units)
        """

        # ignore masses, these are all the same for disk particles
        # m = self.data['m']
        pos, v = self.center_com(com_p, com_v)
        # p = m * v # linear momentum

        if r_lim is not None:
            r = norm(pos, axis=0)
            central = np.where(r < r_lim)
            pos = (pos.T[central]).T
            v = (v.T[central]).T

        # angular momentum of each particle; use first dimension of each array
        L_i = np.cross(pos, v, 0, 0) 

        return np.sum(L_i, axis=0), pos, v

    def rotate_frame(self, to_axis=None, com_p=None, com_v=None, r_lim=None):
        """
        Arg: to_axis (3-vector)
                Angular momentum vector will be aligned to this (default z-hat)

        Returns: (positions, velocities), two arrays of shape (3, N)
                New values for every particle. `self.data` remains unchanged.
        
        Based on Rodrigues' rotation formula
        Ref: https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        """

        if to_axis is None:
            to_axis = np.array([0, 0, 1])
        else:
            to_axis /= norm(to_axis) # we need a unit vector

        L, pos, v = self.angular_momentum(com_p, com_v, r_lim)
        L /= norm(L) # we need a unit vector

        # cross product between L and new axis
        k_vec = np.cross(L, to_axis) # 3-vector
        s_sq = np.sum(k_vec**2) # scalar, sin^2 theta

        # dot product between L and new axis 
        c = np.dot(L, to_axis) # scalar, cos theta

        # rotation matrix, 3x3
        kx, ky, kz = k_vec
        K = np.array([[0, -kz, ky], [kz, 0, -kx], [-ky, kx, 0]])
        R = np.eye(3) + K + K@K * (1 - c) / s_sq

        # Rotate coordinate system
        pos = np.dot(R, pos)
        vel = np.dot(R, v)

        return pos, vel


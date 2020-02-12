# import modules
import numpy as np
from numpy.linalg import norm

import astropy.units as u
import astropy.table as tbl

from galaxy import Galaxy

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
    
    
    def __init__(self, gal, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # subset with our particle type
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
    
    
    def com_p(self, delta=0.1):
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
        # will re-start at half that radius (reduced radius)                                                           
        r_max = max(r_new)/2.0
        
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
            r_max = r_max/2.0
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

        Args: X, Y, Z positions of the COM

        Returns: 3-Vector of COM velocities
        """
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        rv_max = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        xyz_v = self.xyz * u.kpc - xyz_com[:, np.newaxis]
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
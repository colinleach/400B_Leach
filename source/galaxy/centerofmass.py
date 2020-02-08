# import modules
import numpy as np
from numpy.linalg import norm

import astropy.units as u
import astropy.table as tbl

# the next bit is a bodge,
# it may be better to install the package with setup.py
try:
    from galaxy import Galaxy  # works locally
except BaseException:
    from .galaxy import Galaxy  # works on RTD

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


    def COMdefine(self, xyz, m):
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
    
    
    def COM_P(self, delta):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance)                                                                                         
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        xyzCOM = self.COMdefine(self.xyz, self.m)
        # compute the magnitude of the COM position vector.
        RCOM = norm(xyzCOM)

       # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        xyzNew = self.xyz - xyzCOM[:, np.newaxis]
        RNEW = norm(xyzNew, axis=0)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/2.0
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
            # select all particles within the reduced radius 
            # (starting from original x,y,z, m)
            index2 = np.nonzero(RNEW < RMAX)[0]
            xyz2 = self.xyz[:, index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            xyzCOM2 = self.COMdefine(xyz2, m2)
            # compute the new 3D COM position
            RCOM2 = norm(xyzCOM2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            # print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/2.0
            # check this.                                                                                              
            # print ("RMAX", RMAX)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            xyzNew = self.xyz - xyzCOM2[:, np.newaxis]
            RNEW = norm(xyzNew, axis=0)

            # set the center of mass positions to the refined values                                                   
            xyzCOM = xyzCOM2
            RCOM = RCOM2

        # set the correct units using astropy and round all values
        # and then return the COM position vector
        return np.round(xyzCOM, 2) * u.kpc
        
    
    def COM_V(self, xyzCOM):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        xyzV = self.xyz * u.kpc - xyzCOM[:, np.newaxis]
        RV = norm(xyzV, axis=0)
        
        # determine the index for those particles within the max radius
        indexV = np.where(RV < RVMAX)[0]

        # determine the velocity and mass of those particles within the max radius
        vxyznew = self.vxyz[:, indexV]
        mnew = self.m[indexV]  
        
        # compute the center of mass velocity using those particles
        vxyzCOM = self.COMdefine(vxyznew, mnew)

        # return the COM vector                                                                                        
        return np.round(vxyzCOM, 2) * u.km / u.s
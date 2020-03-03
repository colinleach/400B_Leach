

# Homework 5 Solutions
# Mass Profiles
# G. Besla 



# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass


# Gravitational Constant
# converting G to units of kpc*km^2/s^2/Msun
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) # 4.498768e-6*u.kpc**3/u.Gyr**2/u.Msun


class MassProfile:
    # Class to define the Mass and Rotation Curve of a Galaxy
    
    def __init__(self, galaxy, snap):
    # Initialize the instance of this Class with the following properties:
    # galaxy :  string, e.g. "MW"
    #  snap :  integer, e.g  1
    
        # Determine Filename
        # add a string of the filenumber to the value "000"
        ilbl = '000' + str(snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        # create filenames
        self.filename='%s_'%(galaxy) + ilbl + '.txt'
        
        # read in the file                                                                                             
        self.time, self.total, self.data = Read(self.filename)

        # store the mass, positions, velocities of all particles                                
        self.m = self.data['m']#*u.Msun
        self.x = self.data['x']*u.kpc
        self.y = self.data['y']*u.kpc
        self.z = self.data['z']*u.kpc
    
        # store galaxy name
        self.gname = galaxy
    
    
    
    def MassEnclosed(self, ptype, R):
    # Function that determines the MassEnclosed of particles of a given type
    # input: ptype  Particle type,  1=Halo, 2=Disk, 3=Bulge
    #        R   An Array of Radii within which to compute the mass enclosed. 
    # return: an array with the Mass enclosed (units of Msun)
        
        # Determine the COM position using Disk Particles
        # Disk Particles afford the best centroiding.
        # Create a COM object
        COM = CenterOfMass(self.filename,2)
        # Store the COM position of the galaxy
        # Set Delta = whatever you determined to be a good value in Homework 4.
        GalCOMP = COM.COM_P(0.1)
            
        # create an array to store indexes of particles of desired Ptype                                                
        index = np.where(self.data['type'] == ptype)

        # Store positions of particles of given ptype from the COMP. 
        xG = self.x[index] - GalCOMP[0]
        yG = self.y[index] - GalCOMP[1]
        zG = self.z[index] - GalCOMP[2]
            
        # Compute the mag. of the 3D radius
        rG = np.sqrt(xG**2 + yG**2 + zG**2)
            
        # store mass of particles of a given ptype
        mG = self.m[index]
            
        # Array to store enclosed mass as a function of the input radius array
        Menc = np.zeros(np.size(R))
    
        # set up a while loop that continues until the end of the input radius array
        for i in range(np.size(R)):
            # Only want particles within the given radius
            indexR = np.where(rG < R[i]*u.kpc)
            Menc[i] = np.sum(mG[indexR])*1e10         
        
        # return the array of enclosed mass with appropriate units
        return Menc*u.Msun
        
    
    def MassEnclosedTotal(self, R):
    # Determine the total mass of each galaxy point
        # Input:   R   Radius 
        # Returns: Mass in units Msun. 
            
        # Sum up all the mass of each component.
        Menc = self.MassEnclosed(1,R) + self.MassEnclosed(2,R) + self.MassEnclosed(3,R)
    
        # Recall that M33 only has 2 components!  No bulge
        if (self.gname == 'M33'):
            Menc = self.MassEnclosed(1,R) + self.MassEnclosed(2,R)  
          
        return Menc
    
        
        
    def HernquistMass(self, R, scale, Mhalo):
    # Determine the mass enclosed using Hernquist 1990 Mass profile 
    # Input:   R   Radius  
    #         scale   Scale Length  
    #         Mhalo  Total Halo Mass (Msun)
    # Returns: Mass in units Msun. 
    
        # Hernquist 1990 Mass profile
        return Mhalo*R**2/(R+scale)**2*u.Msun
       
        
        
        
    def CircularVelocity(self, ptype, R):
    # Function that determines the circular speed at a given radius using the mass enclosed
    # Input:  ptype  Particle type
    #         R    ARRAY of Radii 
    # Returns:  An ARRAY of Circular speeds in units of km/s
    # NOTE: this method assumes spherical symmetry - but the disk is not spherically symmetric.
    
        # call MassEnclosed to determine array of masses within R.
        Menc = self.MassEnclosed(ptype,R)
        
        # Determine the circular speed as a function of input radius assuming spherical symmetry 
        # note that rad needs units 
        # This will return units of kpc/Gyr
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)
        
        # Return array Vcirc
        return  Vcirc 
    
    
    def CircularVelocityTotal(self, R):
    # Function that determines the combined circular speed, using all galactic components
    # Input:  R Array of Radii
    # Returns: An Array of Circular speends in units of km/s
    
        # compute the total mass enclosed
        Menc = self.MassEnclosedTotal(R)
        
        # Determine the circular speed as a function of input radius assuming spherical symmetry 
        # note that rad needs units 
        # This will return units of kpc/Gyr
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)
        
        # Return array Vcirc
        return  Vcirc 
    
    
    
    def HernquistVCirc(self, R, scale, Mhalo):
    # Function that determines the theoretical circular speed of a Hernquist 1990 Halo
    # input:  R  radius
    #        scale  Scale Length  
    #         Mhalo  Total Halo Mass (Msun)
    # Returns: Vircular Velocity in units km/s.        
             
        # Store the enclosed mass 
        Menc = self.HernquistMass(R,scale,Mhalo)
    
        # Circular speed using enclosed mass
        # rounded to 2 decimal places 
        # You could also write this analytically without first calling HernquistMass
        Vcirc = np.round(np.sqrt(G*Menc/R/u.kpc),2)
        
        # Returns array Vcirc 
        return Vcirc 
    
    
        
        




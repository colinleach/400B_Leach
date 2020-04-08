

# Homework 4 Solutions
# Center ofMass Position and Velocity
# G. Besla,  E. Patel & R. Li

# Homework 6:  COM_P modified

# In[1]:

# import modules
import numpy as np
# astropy provides unit system for astronomical calculations
import astropy.units as u
# astropy also offer a module for showinng results in a table
import astropy.table as tbl
# import previous HW functions
from ReadFile import Read



class CenterOfMass:
    """ Hold the COM position & velocity of a galaxy at a given snapshot """

    
    
    def __init__(self, filename, ptype):
        """ Initialize this class with relevant data """
    
        # read in the file                                                                                             
        self.time, self.total, self.data = Read(filename)
        #print(self.time)                                                                                              

        #create an array to store indexes of particles of desired Ptype                                                
        self.index = np.where(self.data['type'] == ptype)

        # store the mass, positions, velocities of only the particles of the given type                                
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]



    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically                                          
    # input: array (a,b,c) of positions or velocities and the mass                                                     
    # returns: 3 floats  (the center of mass coordinates)                                                              

        # note: since all particles have the same                                                                      
        # mass, when we consider only one type,                                                                             
        # the below is equivalently np.sum(x)/len(x)                                                                   

        # xcomponent Center of mass                                                                                    
        Acom = np.sum(a*m)/np.sum(m)
        # ycomponent Center of mass                                                                                    
        Bcom = np.sum(b*m)/np.sum(m)
        # zcomponent                                                                                                   
        Ccom = np.sum(c*m)/np.sum(m)
        return Acom, Bcom, Ccom
    
    # Modified For Homework 6
    def COM_P(self, delta, VolDec):
        """ Iteratively determine the COM position """
    # input:                                                                                                           
    #        delta (tolerance)  
    #        VolDec (value with which to decrease RMAX)
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

    
        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine                                                   
        XCOM, YCOM, ZCOM = self.COMdefine(self.x,self.y,self.z,self.m)
        # compute the magnitude of the COM position vector.                                                            
        RCOM = np.sqrt(XCOM**2 + YCOM**2 + ZCOM**2)
        # print('init R', RCOM)                                                                                        


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position                                                                          
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        RNEW = np.sqrt(xNew**2.0 + yNew**2.0 +zNew**2.0)

        # find the max 3D distance of all particles from the guessed COM
        # will re-start at a reduced radius specified by input VolDec                                                          
        RMAX = max(RNEW)/VolDec
        
        # pick an initial estimate for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume.                                
        CHANGE = 1000.0

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta):
          # select all particles within the reduced radius (starting from original x,y,z, m)                         
            index2 = np.where(RNEW < RMAX)
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]

            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius                                                                      
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2,y2,z2,m2)
            # compute the new 3D COM position                                                                          
            RCOM2 = np.sqrt(XCOM2**2 + YCOM2**2 + ZCOM2**2)

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # check this                                                                                               
            # print ("DIFF", diff)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by the specified decrement again                                                                
            RMAX = RMAX/VolDec
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

          # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM                                                                                     
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)



            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                
                                                                                                 
            COMP = [XCOM,YCOM,ZCOM]

        # return the COM positon vector
        # set the correct units using astropy                                                                      
        # round all values 
        return np.round(COMP,2)*u.kpc
    
    

    def COM_V(self, COMX,COMY,COMZ):
        """ Return the COM velocity based on the COM position """
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position                              
        xV = self.x[:]*u.kpc - COMX
        yV = self.y[:]*u.kpc - COMY
        zV = self.z[:]*u.kpc - COMZ
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius                                                
        indexV = np.where(RV < RVMAX)

        # determine the velocity and mass of those particles within the mas radius                                     
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]

        # compute the center of mass velocity using those particles                                                    
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew,vynew,vznew, mnew)

        # create a vector to store the COM velocity                                                                    
        # set the correct units usint astropy                                                                          
        # round all values                                                                                             
        COMV = [VXCOM,VYCOM,VZCOM] 

        # return the COM vector  
        # set the correct units using astropy                                                                          
        # round all values   
        return np.round(COMV,2)*u.km/u.s



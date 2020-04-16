#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib


# # Part A :  Mass to Light Ratios 
# 
# Wolf et al. 2010 
# 
# $M(<R_{half}) = \frac {4}{G}\sigma^2 R_e$
# 
# Where $R_{half}$ = 3D half mass radius 
# and $R_e$ is the 2D half mass radius of stars (observed)
# 
# Determine which of the following two systems are galaxies:
# 
# The system 47 Tuc is observed with:  $\sigma = 17.3$ km/s, $R_e = 0.5$ pc, $L_v \sim 10^5 L_\odot$ 
# 
# The system Willman I is observed with: $\sigma = 4.3$ km/s, $R_e = 25$ pc, $L_v = 10^3 L_\odot$

# In[2]:


G = 4.498768e-6 # units of kpc^3/Gyr^2/Msun


# In[3]:


def WolfMass(sigma, Re):
    """ Wolf mass estimator from Wolf+ 2010
    Input sigma = 1D line of sight velocity dispersion in km/s
           Re = 2D radius enclosing half the stellar mass in pc
    Returns estimate of the dynamical mass within the half light radius in Msun"""
    return 4/G*sigma**2*Re/1000


# In[4]:


# 47 Tuc the globular cluster
# Luminosity  L = 1e5 Lsun
# sigma = 17.3 km/s 
# r_c = 0.5 pc

TucMass = WolfMass(17.3,0.5)
print(f" Mass to Light Ratio 47 Tuc: {np.around(TucMass/1e5,3)}")
# this is globular cluster


# In[5]:


# Willman I 
# Luminsoity Lv = 1e3 Lsun
# Sigma = 4.3 km/s
# r_c  =  25 pc 

WillmanMass = WolfMass(4.3,25)
print(f" Mass to Light Ratio Willman 1: {np.around(WillmanMass/1e3,3):}")
# this is a galaxy


# # Part B :  Stellar to Halo Mass Relation
# 
# Following the work of [Moster et al. 2013 (MNRAS, 428, 3121)](https://ui.adsabs.harvard.edu/abs/2013MNRAS.428.3121M/abstract)
# 
# 
# `Equation 2:`                  $ \frac{m}{M} = 2N \left [ \left ( \frac{M}{M_1} \right)^{-\beta} + \left (\frac{M}{M_1} \right)^{\gamma} \right]$ 
# 
# $m$ = stellar mass, $M$ = halo mass
# 
# `Equation 11:`        log $M_1(z) = M_{10} + M_{11} \frac{z}{z+1} $ 
# 
# `Equation 12:`        $N(z) = N_{10} + N_{11} \frac{z}{z+1} $
# 
# `Equation 13:`         $\beta(z) = \beta_{10} + \beta_{11} \frac{z}{z+1} $
# 
# `Equation 14:`         $\gamma(z) = \gamma_{10} + \gamma_{11} \frac{z}{z+1} $

# # Q1 
# 
# Modify the class below by adding a function that takes the `SHMratio` and returns the stellar mass.

# In[6]:


class AbundanceMatching:
    
    def __init__(self, M, z):
        " input: Halo mass (Msun) and Redshift"
        
        #initializing the parameters:
        self.M = M # Halo Mass in Msun
        self.z = z  # Redshift
        
        
    def logM1(self):
        """eq. 11 of Moster 2013
        input : redshift
        output: M1, characteristic mass in log(Msun)
        """
        M10      = 11.59
        M11      = 1.195 
        return M10 + M11*(self.z/(1+self.z))  
    
    
    def N(self):
        """eq. 12 of Moster 2013
        input: redshift
        output: Normalization for eq. 2
        """
        N10      = 0.0351
        N11      = -0.0247
    
        return N10 + N11*(self.z/(1+self.z))
    
    
    def Beta(self):
        """eq. 13 of Moster 2013
        input: redshift
        output: power of the low mass slope"""
        beta10      = 1.376
        beta11      = -0.826
    
        return beta10 + beta11*(self.z/(1+self.z))
    
    def Gamma(self):
        """eq. 14 of Moster 2013
        input: redshift
        output: power of the high mass slope """
        gamma10      = 0.608
        gamma11      = 0.329
    
        return gamma10 + gamma11*(self.z/(1+self.z))
    
    
    def SHMratio(self):
        """ 
        eq. 2 of Moster + 2013
        Inputs: halo mass M in solar masses (NOT in logspce)
           redshift
        Outputs: Stellar mass to halo mass ratio
        """
        M1 = 10**self.logM1() # Converting characteristic mass to Msun from Log(Msun)
        A = (self.M/M1)**(-self.Beta())  # Low mass end
        B = (self.M/M1)**(self.Gamma())   # High mass end
        Norm = 2*self.N() # Normalization
    
        SHMratio = Norm*(A+B)**(-1)
    
        return SHMratio
    
    
 # Q1: add a function to the class that takes the SHM ratio and returns 
# The stellar mass 

    def StellarMass(self):
        """ using eq. 2 of Moster + 2013 to return stellar mass
        Inputs:  halo mass M in solar masses (NOT in logspce)
            redshift
        Outputs: Stellar mass in Msun"""
    
        return self.M*self.SHMratio()


# # Part C : Plot the Moster Relation
# 
# Reproduce the below figure from Moster + 2013 
# Plot this for z=0, 0.5, 1, 2
# 
# ![mos](./MosterFig.png)

# In[7]:


Mh = np.logspace(10,15,1000) # Logarithmically spaced array


# In[8]:


# Define Instances of the Class for each redshift
MosterZ0 = AbundanceMatching(Mh,0)
MosterZ0_5 = AbundanceMatching(Mh,0.5)
MosterZ1 = AbundanceMatching(Mh,1)
MosterZ2 = AbundanceMatching(Mh,2)


# In[9]:



fig,ax = plt.subplots(figsize=(10,8))


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# Plot z = 0
plt.plot(np.log10(Mh), np.log10(MosterZ0.StellarMass()),linewidth = 5, label='z=0')

plt.plot(np.log10(Mh), np.log10(MosterZ0_5.StellarMass()),linewidth = 5, linestyle='-.', label = 'z=0.5')
plt.plot(np.log10(Mh), np.log10(MosterZ1.StellarMass()),linewidth = 5, linestyle='--', label = 'z=1')
plt.plot(np.log10(Mh), np.log10(MosterZ2.StellarMass()),linewidth = 5, linestyle=':', label='z=2')


# Axes labels 
plt.xlabel('log (M$_h$/M$_\odot$)',fontsize=22) 
plt.ylabel('log (m$_\star$/M$_\odot$)', fontsize=22)

# Legend
plt.legend(loc='lower right',fontsize='x-large')

plt.show()

# At every halo mass, the stellar mass is lower backwards in time.
# The characteristic stellar mass is higher at higher redshift -- so massive things do form early, then the smaller
# systems will catch up 
# Lowest stellar mass is 1e8 


# # Part D
# 
# # Q1
# 
# In traditional models of the Magellanic Clouds (prior to 2010), the LMC is thought to have a halo mass of order $3 \times 10^{10}$ M$_\odot$.  According to LCDM theory, what should be the stellar mass of such a halo?  
# 
# How does this compare against the actual observed stellar mass of the LMC at the present day of $3 \times 10^9$ M$_\odot$ ? 
# 
# What is the $\Lambda$CDM expected halo mass? What is the origin of any discrepancy? 

# In[10]:


# Create a class object
L1 = AbundanceMatching(3e10,0)

# Use LMC object to determine its stellar mass
print(f"Expected Stellar Mass of LMC (1e9 Msun) {np.round(L1.StellarMass()/1e9,3)}")

# Compare against the real stellar mass 
LMCstar = 3e9
print(f"Percentage of observed stellar mass: {np.round(L1.StellarMass()/LMCstar,3)*100}%")


# In[11]:


# So what should the LMC's halo mass be to give a stellar mass of 3e9? 
# Start with a guess
Lhalo = 1.65e11

L2 = AbundanceMatching(Lhalo,0) 
# StellarMass(Mhalo, 0)
print(f"Expected Stellar Mass of LMC (1e9 Msun) {np.round(L2.StellarMass()/1e9,3)}")

# LEsson: the Halo mass relation only returns the halo mass of CENTRALS 


# # Q2
# 
# What is the expected stellar mass of an L* galaxy at z=0? 
# 
# What is the expected stellar mass of an L* galaxy at z = 2? 

# In[12]:


print(f'Log M1, characteristic halo mass at z=0: {MosterZ0.logM1()}')
MstarZ0 = AbundanceMatching(10**MosterZ0.logM1(),0)

print(f'Stellar mass of L* galaxy at z=0 (1e10 Msun) : {np.around(MstarZ0.StellarMass()/1e10,2)}')


# In[13]:


print(f'Log M1, characteristic halo mass at z=2: {MosterZ2.logM1()}')
MstarZ2 = AbundanceMatching(10**MosterZ2.logM1(),0)

print(f'Stellar mass of L* galaxy at z=2 (1e10 Msun) : {np.around(MstarZ2.StellarMass()/1e10,2)}')


# In[ ]:





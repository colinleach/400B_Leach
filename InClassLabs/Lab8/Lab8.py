

# # Lab 8 : Star Formation 



import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt


# # Part A
# 
# Create a function that returns the SFR for a given luminosity (NUV, FUV, TIR, Halpha)
# 
# $Log( {\rm SFR} (M_\odot/year)) = Log(Lx (erg/s)) - Log(Cx)$ 
# 
# Including corrections for dust absorption 
# 
# Kennicutt & Evans 2012 ARA&A Equation 12 and Table 1, 2




# Let's try to reproduce SFRs derived for galaxies from UV luminosities measured with Galex. 
# 
# Using Table 1 from Lee et al. 2009
# https://ui.adsabs.harvard.edu/abs/2009ApJ...706..599L/abstract
# 
# We will use galaxy properties from NED:
# https://ned.ipac.caltech.edu/




#  WLM Dwarf Irregular Galaxy


# In[ ]:


#  N24 Sc galaxy


# # Part B Star formation main sequence
# 
# Write a function that returns the average SFR of a galaxy at a given redshift. 
# 
# What is the average SFR of a MW mass galaxy today? at z=1?
# 
# Plot the SFR main sequence for a few different redshifts.
# 
# 
# From Whitaker 2012:
# 
# log(SFR) = $\alpha(z)({\rm log}M_\ast - 10.5) + \beta(z)$
# 
# $\alpha(z) = 0.7 - 0.13z$
# 
# $\beta(z) = 0.38 + 1.14z - 0.19z^2$

# In[ ]:





# In[ ]:


# MW at z=0


# In[ ]:


# MW at z = 1


# In[ ]:


# create an array of stellar masses


# In[ ]:



fig = plt.figure(figsize=(8,8), dpi=500)
ax = plt.subplot(111)

# add log log plots


# Add axis labels
plt.xlabel('Log (Mstar (M$_\odot$))', fontsize=12)
plt.ylabel('Log(SFR (M$_\odot$/year))', fontsize=12)


#adjust tick label font size
label_size = 12
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')

plt.show()

# # Part C  Starbursts
# 
# What are the star formation rates for :
# 
# LIRG: $10^{11}$ L$_\odot$
# 
# ULIRG: $10^{12} $ L$_\odot$
# 
# HLIRG: $10^{13} $ L$_\odot$

# In[ ]:


# normal galaxies 


# In[ ]:


# LIRGs  


# In[ ]:


# ULIRGs


# In[ ]:


# HLIRGs


# In[ ]:





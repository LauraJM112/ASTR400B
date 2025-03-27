
# # Lab 8 : Star Formation 




import numpy as np
from astropy import units as u
from astropy import constants as const

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm



# # Part A
# 
# Create a function that returns the SFR for a given luminosity (NUV, FUV, TIR, Halpha)
# 
# $Log( {\rm SFR} (M_\odot/year)) = Log(Lx (erg/s)) - Log(Cx)$ 
# 
# Including corrections for dust absorption 
# 
# Kennicutt & Evans 2012 ARA&A Equation 12 and Table 1, 2

def StarFormationRate(L, Type, TIR=0):
    '''
    A function that comuptes the star formation rate of a galazy following
    Kennicutt & Evans 2012 Eq 12 (A&A5 0)

    Parameters
    ----------
    L : float
        luminosity of the galaxy (erg/s)
    Type : string
        the wavelength: "FUV", "NUV", "TIR", "Halpha"
    TIR : float, optional
        the total infrared luminosity (erg/s). The default is 0.

    Returns
    -------
    SFR: float
        Log of the star formation rate (Msun/year)

    '''

    if (Type == "FUV"):
        #pull from tables
        logCx = 43.35 #calibration for LFUV to SFR
        TIRc = 0.46 #correction for dust absorbtion
        
    elif (Type == "NUV"):
        logCx = 41.37
        TIRc = 0.27
        
    elif(Type == "Halpha"):
        logCx = 41.27
        TIRc = 0.0024
    elif(Type == "TIR"):
        logCx = 43.41
        TIRc = 0
        
    else:
        print("WARNING: Missing Wavelength. Enter one of the following: FUV, NUV, Halpha, TIR")
        

    #Correct luminosity for dust using TIRc
    Lcorr = L + TIRc*TIR
    
    #star formation rate
    SFR = np.log10(Lcorr - logCx)
    
    return SFR

# Let's try to reproduce SFRs derived for the WLM Dwarf Irregular Galaxy using 
#UV luminosities measured with Galex. 
# 
# Compare results to Table 1 from Lee et al. 2009 (who used the older Kennicutt 98 methods)
# https://ui.adsabs.harvard.edu/abs/2009ApJ...706..599L/abstract
# 
# We will use galaxy properties from NED (Photometry and SED):
# https://ned.ipac.caltech.edu/



# First need the Luminosity of the Sun in the right units
LsunErgS = const.L_sun.to(u.erg/u.s).value



#  WLM Dwarf Irregular Galaxy
NUV_WLM = 1.71e7*LsunErgS #from NED GALEX data
TIR_WLM = 2.48e6*LsunErgS + 3.21e5*LsunErgS + 2.49e6*LsunErgS
#TIR = NIR + MIR + FIR

#test code
print(StarFormationRate(NUV_WLM, 'NUV', TIR_WLM)) #this is wrong should be -2 ish
#print(StarFormationRate(1e6*LsunErgS, 'blah', 5e6*LsunErgS))
# # Part B Star formation main sequence
# 
# 1) Write a function that returns the average SFR of a galaxy at a given redshift, given its stellar mass
# 
# 2) What is the average SFR of a MW mass galaxy today? at z=1?
# 
# 3) Plot the SFR main sequence for a few different redshifts from 1e9 to 1e12 Msun.
# 
# 
# From Whitaker 2012:
# 
# log(SFR) = $\alpha(z)({\rm log}M_\ast - 10.5) + \beta(z)$
# 
# $\alpha(z) = 0.7 - 0.13z$
# 
# $\beta(z) = 0.38 + 1.14z - 0.19z^2$

# # Step 1
def SFRMainSequence(Mstar, z):
    '''
    A function that computes the average star formation rate of a galaxy as a 
    function of stellar mass and redshift

    Parameters
    ----------
    Mstar : float
        Stellar mass of the galaxy (Msun)
    z : float
        Redshift

    Returns
    -------
    SFR: float
        log of the SFR (Msun/year)

    '''
    alpha = 0.7 - 0.13*z
    beta = 0.38 + 1.14*z - 0.19*z**2
    
    SFR = alpha*(np.log10(Mstar) - 10.5) + beta
    
    return SFR




# # Step 2
# MW at z=0
MWmstar = 7.5e10 #stellar disk mass from simulation
print(SFRMainSequence(MWmstar, 0))
print(10**SFRMainSequence(MWmstar, 0))


# MW at z = 1
print(SFRMainSequence(MWmstar, 1))
print(10**SFRMainSequence(MWmstar, 1))

# # Step 3
# create an array of stellar masses
Mass = np.linspace(1e8, 1e12)




fig = plt.figure(figsize=(8,8), dpi=500)
ax = plt.subplot(111)

# add log log plots
plt.plot(np.log10(Mass), SFRMainSequence(Mass, 0), color='blue', label = 'z=0')
plt.plot(np.log10(Mass), SFRMainSequence(Mass, 1), color='red', label = 'z=1', linestyle = "-.")
plt.plot(np.log10(Mass), SFRMainSequence(Mass, 2), color='green', label = 'z=2', linestyle = ':')
plt.plot(np.log10(Mass), SFRMainSequence(Mass, 3), color='black', label = 'z=3', linestyle = "--")
# Add axis labels
plt.xlabel('Log(Mstar (M$_\odot$))', fontsize=12)
plt.ylabel('Log(SFR (M$_\odot$/year))', fontsize=12)


#adjust tick label font size
label_size = 12
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size

# add a legend with some customizations.
legend = ax.legend(loc='upper left',fontsize='x-large')
plt.show()

# Save file
#plt.savefig('Lab8_SFR_MainSequence.png')


# # Part C  Starbursts
# 
# Use your `StarFormationRate` code to determine the typical star formation rates for the following systems with the listed Total Infrared Luminosities (TIR): 
# 
# Normal Galaxies: $10^{10}$ L$_\odot$
# 
# LIRG: $10^{11}$ L$_\odot$
# 
# ULIRG: $10^{12} $ L$_\odot$
# 
# HLIRG: $10^{13} $ L$_\odot$



# normal galaxies 
TIR_Normal = 1e10*LsunErgS
print(10**StarFormationRate(TIR_Normal, "TIR"))



# LIRGs  




# ULIRGs




# HLIRGs

TIR_HLIRG = 1e13*LsunErgS
print(10**StarFormationRate(TIR_HLIRG, "TIR"))





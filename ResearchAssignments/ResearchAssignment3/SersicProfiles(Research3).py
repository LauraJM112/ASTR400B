#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 30 13:40:12 2025

@author: lauramack
"""

## Research Assignment 3
#How well is the remnant described as a classical elliptical galaxy based on 
#its surface density profile? What is the best fit sersic profile? 

# Find the Sersic profile for the merged galaxy to find its morphology

'''
Import sersic lab, data from post-merger snapshot
find COM and correct positions
change to cylindrical (like in lab6)
find surface density, use function from lab 6
Make function for mean squared error
Initalize array to store index, profile and error
Set up loop over indices (0.5-10?)
Call MSE
Store in array
Save array to file
Plot profiles for all indices
Find smallest error and plot alongside data
'''

##import modules
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
#from Lab6 import sersicE
#from Lab6 import SurfaceDensity
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile
from GalaxyMass import ComponentMass


#find COM and correct positions
M31COM = CenterOfMass("M31_000.txt", 3) #bulge particles
# Center of Mass of M31 
M31_COM_p = M31COM.COM_P(0.1)

#correct to center of mass and store x,y,z,m
x = M31COM.x - M31_COM_p[0].value
y = M31COM.y - M31_COM_p[1].value
z = M31COM.z - M31_COM_p[2].value
m = M31COM.m # units of 1e10

#change to cylindrical (like in lab6)
cyl_r = np.sqrt(x**2 + y**2) # radial
cyl_theta = np.arctan2(y,x) # theta


#find surface density (this is all pulled from lab 6)
def SurfaceDensity(r,m):
    """ Function that computes the surface mass density profile
    given an array of particle masses and radii 
     
    PARMETERS
    ---------
        r : array of `floats` - cyclindrical radius [kpc]
        m : array of `floats` - particle masses [1e10 Msun] 
    
    RETURNS
    -------
        r_annuli : array of `floats` -  radial bins for the 
            annuli that correspond to the surface mass density profile
            
        sigma: array of `floats` - surface mass density profile 
         [1e10 Msun/kpc^2] 
        
    
    """
    
    # Create an array of radii that captures the extent of the bulge
    # 95% of max range of bulge
    radii = np.arange(0.1, 0.95 * r.max(), 0.1)

    # create a mask to select particles within each radius
    # np.newaxis creates a virtual axis to make cyl_r_mag 2 dimensional
    # so that all radii can be compared simultaneously
    # a way of avoiding a loop - returns a boolean 
    enc_mask = r[:, np.newaxis] < radii

    # calculate mass of bulge particles within each annulus.  
    # relevant particles will be selected by enc_mask (i.e., *1)
    # outer particles will be ignored (i.e., *0)
    # axis =0 flattens to 1D
    m_enc = np.sum(m[:, np.newaxis] * enc_mask, axis=0)

    # use the difference between m_enc at adjacent radii 
    # to get mass in each annulus
    m_annuli = np.diff(m_enc) # one element less then m_enc
    
    
    # Surface mass density of stars in the annulus
    # mass in annulus / surface area of the annulus. 
    # This is in units of 1e10
    sigma = m_annuli / (np.pi * (radii[1:]**2 - radii[:-1]**2))
    # array starts at 0, but here starting at 1 and
    # subtracting radius that ends one index earlier.
    
    # Define the range of annuli
    # here we choose the geometric mean between adjacent radii
    r_annuli = np.sqrt(radii[1:] * radii[:-1]) 

    return r_annuli, sigma

r_annuli, sigmaM31bulge = SurfaceDensity(cyl_r, m)

#find mass profile
M31mass = MassProfile("M31",700) 
# Determine the disk mass profile
# use the annuli defined for the surface mass density profile
bulge_mass = M31mass.MassEnclosed(3,r_annuli).value
# Determine the total mass of the bulge
bulge_total = ComponentMass("M31_000.txt", 3)*1e12
# Half the total bulge mass
b_half = bulge_total/2.0
# Find the indices where the bulge mass is larger than b_half
index = np.where(bulge_mass > b_half)
# Define the Effective radius of the bulge
re_bulge = r_annuli[index][0]*3/4

#get sersic profile (function pulled from lab 6)
def sersicE(r, re, n, mtot):
    """ Function that computes the Sersic Profile for an Elliptical 
    System, assuming M/L ~ 1. As such, this function is also the 
    mass surface density profile. 
    
    PARMETERS
    ---------
        r: `float`
            Distance from the center of the galaxy (kpc)
        re: `float`
            The Effective radius (2D radius that contains 
            half the light) (kpc)
        n:  `float`
            sersic index
        mtot: `float`
            the total stellar mass (Msun)

    RETURNS
    -------
        I: `array of floats`
            the surface brightness or mass density (M/L = 1)
            profile for an elliptical in Lsun/kpc^2

    """
    # M/L = 1 
    lum = mtot
    
    # the effective surface brightness
    Ie = lum/7.2/np.pi/re**2
    
    # break down the sersic profile
    a = (r/re)**(1/n)
    b = -7.67*(a-1)
    
    # the surface brightness or mass density profile
    I= Ie*np.exp(b) 
    
    return I


SersicM31Bulge = sersicE(r_annuli, re_bulge, 4, bulge_total)


# Make function for mean squared error
#MSE = 1/n * sum(obs - predicted)**2
def MSE(Sersic, sim_data):
    '''
    A function to find the mean squared error between observed (sersic profile)
    data and the predicted (data from simulation)
    uses MSE = 1/n*sum(obs - predicted)**2

    Parameters
    ----------
    Sersic : array of floats
        the sersic profile
    sim_data : array of floats
        the expected surface density from the simulation

    Returns
    -------
    MSE : float
        the error between the two arrays

    '''
    n = len(sim_data) #number of data points
    print("n =",n)
    diff = np.subtract(Sersic, sim_data) #find difference
    diff_square = diff**2 #square it
    MSE = (1/n)* np.sum(diff_square) #take sum and multiply by 1/n
    
    return MSE

fig, ax = plt.subplots(figsize=(9, 8))

#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Surface Density Profile
# YOU ADD HERE
ax.loglog(r_annuli, sigmaM31bulge, lw=2, label='Sim bulge')

# Sersic fit to the surface brightness Sersic fit
# YOU ADD HERE
ax.loglog(r_annuli, SersicM31Bulge/1e10, linestyle='-.', lw=3,
         label = 'Sersic n=4')


plt.xlabel('log r [ kpc]', fontsize=22)

# note the y axis units
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]', 
          fontsize=22)

plt.title('M31 Bulge', fontsize=22)

#set axis limits
plt.xlim(1,50)
plt.ylim(1e-5,0.1)

ax.legend(loc='best', fontsize=22)
fig.tight_layout()
print(MSE(SersicM31Bulge, sigmaM31bulge))
print(len(sigmaM31bulge))
# Initalize array to store index, profile and error
# Set up loop over indices (0.5-10?)
# Call MSE
# Store in array
# Save array to file
# Plot profiles for all indices
# Find smallest error and plot alongside data


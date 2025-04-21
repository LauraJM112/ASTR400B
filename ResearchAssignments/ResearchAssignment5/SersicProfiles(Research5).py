#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 11:45:17 2025

@author: lauramack
"""
#import modules
import numpy as np
import matplotlib.pyplot as plt

#my modules
from COM import COM_P #my version of Center of Mass that just runs COM_define and COM_P
from COM import Position #reads in files and return x,y,z,m which get passed to COM_P
from Lab6func import SurfaceDensity # put the functions from lab 6 in own file
from Lab6func import sersicE
from Lab6func import MassEnclosed

#from GalaxyMass import ComponentMass

#%%
x,y,z,m = Position() #get disk/bulge particles
COM = COM_P(x, y, z, m) #get center of mass (units of kpc)

#correct positions for COM
x = x - COM[0]
y = y - COM[1]
z = z - COM[2]
r = np.sqrt(x**2+y**2+z**2)

#change to cylindrical (like in lab6)
cyl_r = np.sqrt(x**2 + y**2) # radial
cyl_theta = np.arctan2(y,x) # theta

#run surface density
r_annuli, sigma = SurfaceDensity(cyl_r, m)
print(r_annuli)

#%%
elliptical_mass = MassEnclosed(r_annuli, r, m)


# Determine the total mass of remnant
elliptical_total = sum(m)


# Half the total mass
e_half = elliptical_total/2.0
print("half elliptical:", f"{e_half:.2e}")
# Find the indices where the merger mass is larger than e_half
index = np.where(elliptical_mass > e_half)
print("elliptical_mass[index][0]", f"{elliptical_mass[index][0]:.2e}")
# Define the Effective radius of the bulge
re_elliptical = r_annuli[index][0]*3/4
print("re_elliptical", re_elliptical)

#get sersic profile
SersicP = sersicE(r_annuli, re_elliptical, 4, elliptical_total)





#%%
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
    diff = np.subtract(Sersic, sim_data) #find difference
    diff_square = diff**2 #square it
    MSE = (1/n)* np.sum(diff_square) #take sum and multiply by 1/n
    
    return MSE



#%%
### Plot De Vaucoulers profile
plt.loglog(r_annuli, sigma, label = "simulation")
plt.loglog(r_annuli, SersicP, label = "Sersic n=4", ls='-.')
plt.xlim((0.1, 70))
plt.title("Merger Remnant with De Vaucouleurs profile")
plt.xlabel('log r [ kpc]')
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]')
plt.legend()
plt.show()

#%%
#loop over Sersic indices
ind = np.arange(2, 6, 0.01) # get indicies to test
Sersic = np.zeros((len(ind), len(r_annuli))) #set up place to store values
for i in range(0, len(ind)): #loop through indicies and find sersic profiles
    Sersic[i] = sersicE(r_annuli, re_elliptical, ind[i], elliptical_total)
    

#%%
error = np.zeros((len(ind), 2)) #make array to store errors

for i in range(len(Sersic)): #loop over sersic array
    ind_error = MSE(Sersic[i], sigma) #run MSE for each profile
    error[i] = [ind[i], ind_error]#store in error array
res = min(error, key=min) # list with minimum error value and index

#print(error)

print("Minimum element sublist is : " + str(res))

#%%
plt.plot(error[:,0], error[:,1]) #plot all error measurments
plt.title("MSE error as a function of Sersic Index")
plt.scatter(res[0], res[1], color='tab:orange') #plot point at minimum error 
plt.xlabel("Index")
plt.ylabel("Error")
plt.ylim(0, 0.007)
plt.xlim(2.5, 5)
plt.show()

#%%
best_profile = sersicE(r_annuli, re_elliptical, res[0], elliptical_total)
best_index = np.round(res[0], 3)
plt.loglog(r_annuli, sigma, label = "simulation")
plt.loglog(r_annuli, best_profile, label = f"Sersic n={best_index}", ls='-.')
plt.xlim((0.1, 70))
plt.title("Merger Remnant with best fit profile")
plt.xlabel('log r [ kpc]')
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]')
plt.legend()
plt.show()
    

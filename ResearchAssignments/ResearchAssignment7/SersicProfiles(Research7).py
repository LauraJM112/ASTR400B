#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 11:45:17 2025

@author: lauramack
"""
#import modules
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

#my modules
from COM import COM_P #my version of Center of Mass that just runs COM_define and COM_P
from COM import Remnant_Position #reads in files and return x,y,z,m which get passed to COM_P
from COM import MW_Position #reads in MW_000 and gets disk and bulge particles
from Lab6func import SurfaceDensity # put the functions from lab 6 in own file
from Lab6func import sersicE 
from Lab6func import MassEnclosed 
from Lab6func import ChiSquared # Chi squared error, written specifically for this

#%%
x,y,z,m = MW_Position() #Either MW_Position() or Remnant_Positions() depending on which state is of interest
COM = COM_P(x, y, z, m) #get center of mass (units of kpc)

#correct positions for COM
x = x - COM[0]
y = y - COM[1]
z = z - COM[2]
r = np.sqrt(x**2+y**2+z**2)

#change to cylindrical 
cyl_r = np.sqrt(x**2 + y**2) # radial
cyl_theta = np.arctan2(y,x) # theta

#run surface density
r_annuli, sigma = SurfaceDensity(cyl_r, m)
#print(r_annuli)

#%%
elliptical_mass = MassEnclosed(r_annuli, r, m)

# Determine the total mass 
elliptical_total = sum(m)


# Half the total mass
e_half = elliptical_total/2.0
#print("half elliptical:", f"{e_half:.2e}")
# Find the indices where the merger mass is larger than e_half
index = np.where(elliptical_mass > e_half)
#print("elliptical_mass[index][0]", f"{elliptical_mass[index][0]:.2e}")
# Define the Effective radius of the bulge
re_elliptical = r_annuli[index][0]*3/4
#print("re_elliptical", re_elliptical)

#get sersic profile for n=4
Sersic4 = sersicE(r_annuli, re_elliptical, 4, elliptical_total)

#%%
#loop over Sersic indices
ind = np.arange(1, 10, 0.01) # get indicies to test
Sersic = np.zeros((len(ind), len(r_annuli))) #set up place to store values
for i in range(0, len(ind)): #loop through indicies and find sersic profiles
    Sersic[i] = sersicE(r_annuli, re_elliptical, ind[i], elliptical_total)
    

#%%
error = np.zeros((len(ind), 2)) #make array to store errors
min_error = 1e5 #set big number to start
best_index = 0

for i in range(len(Sersic)): #loop over sersic array

    ind_error = ChiSquared(Sersic[i], sigma) #run chi squared for each profile
    error[i] = [ind[i], ind_error]#store in error array
    if ind_error < min_error:
        min_error = ind_error
        best_index = ind[i]
 

res = [best_index, min_error] #get minimum error and index
# printing the result
print ("Minimum element sublist is : " + str(res))


best_index = np.round(res[0], 3)

#%%
########
# Plot things
########
plt.semilogy(error[:,0], error[:,1]) #plot all error measurments

plt.suptitle("Chi Squared error as a function of Sersic Index", size=12)
plt.title(f"Minimum error is at index {best_index} with value of {round(res[1], 2)}", size=9)
plt.scatter(res[0], res[1], color='tab:orange') #plot point at minimum error 
plt.xlabel("Index")
plt.ylabel("Error")
plt.savefig("Plots/MW_error")
plt.show()


#%%
best_profile = sersicE(r_annuli, re_elliptical, best_index, elliptical_total) #best fit
deVaucouleurs = sersicE(r_annuli, re_elliptical, 4, elliptical_total) #classical elliptical fit

plt.loglog(r_annuli, sigma, label = "simulation")
plt.loglog(r_annuli, best_profile, label = f"Best fit profile n={best_index}", ls='-.')
plt.loglog(r_annuli, deVaucouleurs, label = "de Vaucouleurs profile n=4", ls =':')
plt.xlim((0.1, 32)) #32 is good for MW, 70 for remnant
plt.ylim(1e-6, 10)
plt.title("MW with best fit profile")
plt.xlabel('log r [ kpc]')
plt.ylabel(r'log $\Sigma_{bulge}$ [$10^{10} M_\odot$ / kpc$^2$]')
plt.legend()
plt.savefig("Plots/Best profile for MW")
plt.show()



    

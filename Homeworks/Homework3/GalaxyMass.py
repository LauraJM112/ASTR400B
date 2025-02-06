#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb  2 13:50:47 2025

@author: lauramack
"""

from ReadFile import Read
import astropy.units as u
import numpy as np

def ComponentMass(filename, part_type):
    '''
    This function computes the mass of a part of a galaxy. 

    Parameters
    ----------
    filename : string
        Name of file in question
    part_type : integer (1,2,3)
        The area of interest in the simulation
    Returns
    -------
    TotalMass: float
        The mass present in that portion of the galaxy

    '''
    time, num, data = Read(filename) # run Read 
    type_index = np.where(data["type"] == part_type) #index data for particles  
        # in area of interest   
    type_data = data[type_index] # filter data to only match the type
    Masses = type_data["m"] # get mass in units of 1e10  Msun
    MassMSun = Masses *10**10 #convert to units of Msun
    Mass12 = MassMSun/(10**12) #convert to units of 10**12 Msun
    
    TotalMass = np.round(sum(Mass12), 3) #add up masses and round to 3 places
    
    return TotalMass #in units of 1e12 Msun


#these are for computing things and I should delete
filename = 'MW_000.txt'
  

halo = ComponentMass(filename, 1)
disk = ComponentMass(filename, 2)
bulge = ComponentMass(filename, 3)
total = halo + disk + bulge

print(halo, "&", disk, "&", bulge, "&", total)

stellar = np.round(bulge + disk, 3)
print(np.round(stellar/total, 3))



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 11:32:31 2025

@author: lauramack
"""
import astropy.units as u
import numpy as np
from ReadFile import Read

def ParticleInfo(filename, part_type, num_particles):
    '''
    This function 

    Parameters
    ----------
    filename : string
        the name of a file from MW simulation
    part_type : string
        the type of object of inerest (ex: halo, disk)
    num_particles : TYPE
        DESCRIPTION.

    Returns
    -------
    distance : Array
        The distances of each particle in the area of inerest to the viewer 
        in units of kpc
    velocity : Array
        The velocities of each particle in the area of interest in units 
        of km/s
    mass : Array
        The mass of each particle in units of Mass of the sun

    '''
    time, num_particles, data = Read(filename) #get values from read function
    
    type_index = np.where(data["type"] == part_type) #index data for particles  
        # in area of interest
    type_data = data[type_index] # filter data to only match the type
    
    
    x = (type_data["x"]) #grab x values
    y = type_data["y"] #grab y values
    z = type_data["z"] #grab z values
    distance = (x**2 + y**2 + z**2)**0.5*u.kpc #use pythagorean theorem to 
        #find the magnitude of the distance
    distance = np.around(distance, 3) #round to 3 places
    distance = distance.to(u.lightyear)
    
    
    vx = type_data["vx"] #grab x velocities
    vy = type_data["vy"] #grab y velocities
    vz = type_data["vz"] #grab z velocities
    velocity = (vx**2 + vy**2 + vz**2)**0.5*u.km/u.s #use pythagorean theorem to
        #find the magnitude of the velocity
    velocity = np.around(velocity,3) #round to 3 places
    
    mass = type_data["m"]*u.Msun #grab mass values and set to units of M_sun
    
    
    return distance, velocity, mass

    

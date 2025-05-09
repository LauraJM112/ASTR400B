#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 16 11:07:12 2025

@author: lauramack
"""
import numpy as np
from ReadFile import Read

def MW_Position():
    time, num, data = Read("MW_000.txt")
    
    #get inidices for disk and bulge paricles
    disk_index = np.where(data['type'] == 2)
    bulge_index = np.where(data['type'] == 3)
    
    #pull only disk/bulge particles from data
    disk_particles = data[disk_index]
    bulge_particles = data[bulge_index]
    
    total_data = np.concatenate((disk_particles, bulge_particles)) #put disk/bulge into 1 array
    
    x = total_data["x"] #grab x values in kpc
    y = total_data["y"] #grab y values
    z = total_data["z"] #grab z values
    m = total_data["m"] #get masses in 1e10 msun
    
    return x,y,z,m

def Remnant_Position():
    """
    Reads in MW and M31 at snapshot 500 gives their position and mass information
    for the disk and bulge

    Returns
    -------
    x : array
        the x position of disk/bulge particles
    y : array
        the  yposition of disk/bulge particles
    z : array
        the z position of disk/bulge particles
    m : array
        the  mass of the bulge/disk particles

    """
    #read in files
    M31time, M31num, M31data = Read("M31_500.txt")
    MWtime, MWnum, MWdata = Read("MW_500.txt")

    #M31
    #get the indices for disk, bulge particles
    M31_disk_index = np.where(M31data['type'] == 2)
    M31_bulge_index = np.where(M31data['type'] == 3)
    #grab values at those indices
    M31_disk = M31data[M31_disk_index]
    M31_bulge = M31data[M31_bulge_index]

    #MW
    #get the indices for disk, bulge particles
    MW_disk_index = np.where(MWdata['type'] == 2)
    MW_bulge_index = np.where(MWdata['type'] == 3)
    #grab values at those indices
    MW_disk = MWdata[MW_disk_index]
    MW_bulge = MWdata[MW_bulge_index]

    #put all of that together
    data = np.concatenate((M31_disk, M31_bulge, MW_disk, MW_bulge))

    x = data["x"] #grab x values in kpc
    y = data["y"] #grab y values
    z = data["z"] #grab z values
    m = data["m"] #get masses in 1e10 msun
    
    return x,y,z,m

#pulled from center of mass class
def COMdefine(a,b,c,m):
    ''' Method to compute the COM of a generic vector quantity by direct weighted averaging.
    
    PARAMETERS
    ----------
    a : `float or np.ndarray of floats`
        first vector component
    b : `float or np.ndarray of floats`
        second vector component
    c : `float or np.ndarray of floats`
        third vector component
    m : `float or np.ndarray of floats`
        particle masses
    
    RETURNS
    -------
    a_com : `float`
        first component on the COM vector
    b_com : `float`
        second component on the COM vector
    c_com : `float`
        third component on the COM vector
    '''

    # write your own code to compute the generic COM 
    #using Eq. 1 in the homework instructions
    masses = np.sum(m)
    # xcomponent Center of mass
    a_num = np.sum(a*m) #numerator 
  
    a_com = a_num/masses
    # ycomponent Center of mass
    b_num = np.sum(b*m)
    b_com = b_num/masses
    # zcomponent Center of mass
    c_num = np.sum(c*m)
    c_com = c_num/masses
    
    # return the 3 components separately
    return a_com, b_com, c_com

#pulled from COM
def COM_P(x, y, z, m, delta=0.1):
    '''Method to compute the position of the center of mass of the galaxy 
    using the shrinking-sphere method.

    PARAMETERS
    ----------
    delta : `float, optional`
        error tolerance in kpc. Default is 0.1 kpc
    
    RETURNS
    ----------
    p_COM : `np.ndarray of astropy.Quantity'
        3-D position of the center of mass in kpc
    '''                                                                     

    # Center of Mass Position                                                                                      
    ###########################                                                                                    

    # Try a first guess at the COM position by calling COMdefine                                                   
    x_COM, y_COM, z_COM = COMdefine(x, y, z, m)
    # compute the magnitude of the COM position vector.
    # write your own code below
    r_COM = (x_COM**2 + y_COM**2 + z_COM**2)**0.5 # magnitude of COM


    # iterative process to determine the center of mass                                                            

    # change reference frame to COM frame                                                                          
    # compute the difference between particle coordinates                                                          
    # and the first guess at COM position

    x_new = x - x_COM
    y_new = y - y_COM
    z_new = z - z_COM
    r_new = (x_new**2 +y_new**2+z_new**2)**0.5

    # find the max 3D distance of all particles from the guessed COM                                               
    # will re-start at half that radius (reduced radius)                                                           
    r_max = max(r_new)/2.0
    
    # pick an initial value for the change in COM position                                                      
    # between the first guess above and the new one computed from half that volume
    # it should be larger than the input tolerance (delta) initially
    change = 1000.0

    # start iterative process to determine center of mass position                                                 
    # delta is the tolerance for the difference in the old COM and the new one.    
    
    while (change > delta):
        # select all particles within the reduced radius (starting from original x,y,z, m)
        # write your own code below (hints, use np.where)
        
        #r of all particles (I think)
       # r_ind = (self.x**2 + self.y**2 + self.z**2)**0.5 
        
        index2 = np.where(r_new <= r_max) 
        x2 = x[index2]
        y2 = y[index2]
        z2 = z[index2]
        m2 = m[index2]
        
        #r2 = (x2**2 + y2**2 + x2**2)**0.5

        # Refined COM position:                                                                                    
        # compute the center of mass position using                                                                
        # the particles in the reduced radius
        # write your own code below
        x_COM2, y_COM2, z_COM2 = COMdefine(x2, y2, z2, m2)
        # compute the new 3D COM position
        r_COM2 = (x_COM2**2 + y_COM2**2 + z_COM2**2)**0.5

        # determine the difference between the previous center of mass position                                    
        # and the new one.                                                                                         
        change = np.abs(r_COM - r_COM2)
        # uncomment the following line if you want to check this                                                                                               
        #print ("CHANGE = ", change)                                                                                     

        # Before loop continues, reset : r_max, particle separations and COM                                        

        # reduce the volume by a factor of 2 again                                                                 
        r_max /= 2.0
        # check this.                                                                                              
        #print ("maxR", r_max)                                                                                      

        # Change the frame of reference to the newly computed COM.                                                 
        # subtract the new COM
        # write your own code below
        x_new = x - x_COM2
        y_new = y - y_COM2
        z_new = z - z_COM2
        r_new = (x_new**2+y_new**2+z_new**2)**2

        # set the center of mass positions to the refined values                                                   
        x_COM = x_COM2
        y_COM = y_COM2
        z_COM = z_COM2
        r_COM = r_COM2

        # create an array (np.array) to store the COM position                                                                                                                                                       
        p_COM = np.array([x_COM, y_COM, z_COM])
        p_COM = np.round(p_COM, 3) #round
        #p_COM = p_COM*u.kpc #convert to astropy quantity

    # set the correct units using astropy and round all values
    # and then return the COM positon vector
    return p_COM
    





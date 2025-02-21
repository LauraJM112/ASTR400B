#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:46:50 2025

@author: lauramack
"""

## Homework 5

#import modules
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.constants import G
from ReadFile import Read
from CenterOfMass import CenterOfMass


class MassProfile:
    # a class to find the mass profile of the galaxies 
    
    def __init__(self, galaxy, Snap):
        # add a string of the filenumber to the value “000” 
        ilbl = '000' + str(Snap)
        # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename="%s"%(galaxy) +"_"+ ilbl + '.txt'
        
        #read in data
        self.time, self.total, self.data = Read(self.filename)
        
        self.m = self.data['m'] #grab mass

        self.x = (self.data["x"]) #grab x values
        self.y = self.data["y"] #grab y values
        self.z = self.data["z"] #grab z values
        
        #get velocities
        self.vx = self.data["vx"]
        self.vy = self.data["vy"]
        self.vz = self.data["vz"]
        
        from astropy.constants import G
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #put G in right units
        
        #store name of galaxy
        self.gname = galaxy
        
    def MassEnclosed(self, p_type, radii):
        '''
        finds the mass enclosed by a set radius from the center of mass

        Parameters
        ----------
        p_type : integer
            the type of particles in question
            (halo, disk, bulge)
        radii : array of floats
            the various radii in question

        Returns
        -------
        total_mass : array of astropy quantities
            the mass enclosed at each radius

        '''
        type_index = np.where(self.data["type"] == p_type) #index data for particles  
            # in area of interest
        type_data = self.data[type_index] # filter data to only match the type
        #get center of mass of the galaxy
        gal_CM = CenterOfMass(self.filename, 2)
        COM_p = gal_CM.COM_P(0.1)  
        size = len(radii) #size of radii array
        total_mass = np.zeros(size) #initalize array to store the total mass enclosed at each radii
        x = type_data["x"] - COM_p[0].value #adjust positions to be in COM frame
        y = type_data["y"] - COM_p[1].value
        z = type_data["z"] - COM_p[2].value
        r = (x**2 + y**2 + z**2)**0.5*u.kpc #get magnitude of distance for each particle

        r_COM = r #- COM_p #radius from center of mass
        
        for i in range(size): #loop over radii
            index = np.where(r_COM.value <= radii[i]) #get index for each 
            masses = type_data["m"][index] #grab mass of each particle of the specified type enclosed
            mass = sum(masses)
            total_mass[i] = mass #store in total_mass
            
        total_mass = total_mass*1e10*u.Msun #add units of solar masses
        return total_mass
    
    
    
    def TotalMassEnclosed(self, radii):
        '''
        Uses Mass Enclosed fuction to find the total mass enclosed for all 
        particle types at various radii

        Parameters
        ----------
        radii : an array of floats
            The radii where the mass enclosed will be calculated (units of kpc)

        Returns
        -------
        Mtotal : an array of astropy quantities 
            the mass enclosed for all particles at the specified radii (units of solar mass)

        '''
        Mhalo = self.MassEnclosed(1, radii) #halo mass
        Mdisk = self.MassEnclosed(2, radii) #disk mass
        if self.gname == "M33": #M33 does not have a bulge so sets bulge mass to 0
            Mbulge = 0
        else:
            Mbulge = self.MassEnclosed(3, radii) #bulge mass
        
        Mtotal = Mhalo + Mdisk + Mbulge #add masses together
        
        return Mtotal
    
    
    def HernquistMass(self, radius, a, Mhalo):
        '''
        Find the mass enclosed at a set radius using a hernquist profile
        M = Mhalo * r**2/(a+r)**2
        
        Parameters
        ---------
        radius: float 
            The radius where mass enclosed should be calculated
        a: float
            The scale factor
        Mhalo: astropy quanity (Msun) 
            Mass of the halo
            
        Returns
        -------
        HernMass: astropy quantity
            mass enclosed given by herquist profile
        
        '''
        num = Mhalo * radius**2 #numerator
        den = (a+radius)**2 #denominator
        HernMass = num/den * u.Msun #divide to find hernquist mass and sets to units of Msun
        
        return HernMass
    
    
    def CircularVelocity(self, p_type, radii):
        """
        Calculates the circular velocity using sqrt(GM/r**2) 

        Parameters
        ----------
        p_type : integer
            the type of particle of interest
        radii : array of floats 
            the radii where circular velocity should be calculated

        Returns
        -------
        Vcirc : astropy quanity (km/s)
            the circular velocity

        """
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #put G in right units
        M = self.MassEnclosed(p_type, radii) #get enclosed mass at each radius
        Vcirc = np.sqrt(G*M/radii**2)*u.km/u.s #use formula to compute
        
        return Vcirc


    def CircularVelocityTotal(self, radii):
        '''
        Calculates the total circular velocity using total enclosed mass 
        and sqrt(GM/r**2)
        
        Parameters
        ---------
        radii: array of floats
            the radii where circualr velocity is to be computed
            
        Returns
        -------
        VCircTot: astropy quantity (km/s)
            total circular velocity
        '''

        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #put G in right units
        Mass = self.TotalMassEnclosed(radii) #get total enclosed mass
        VCircTot = np.sqrt(G*Mass/radii**2)*u.km/u.s #use formula to compute
        
        return VCircTot
        
    def HernquistVCirc(self, radii, a, Mhalo):
        '''
        calculates the circular velocity using the hernquist profile

        Parameters
        ----------
        radii : array of floats
            radii of interest
        a : float
           scale factor for hernquist profile
        Mhalo : float
            mass of halo for hernquist profile

        Returns
        -------
        VHern : astropy quanity (km/s)
            circular velocity

        '''
        self.G = G.to(u.kpc*u.km**2/u.s**2/u.Msun) #put G in right units
        Mass = self.HernquistMass(radii, a, Mhalo) #get total enclosed mass
        VHern = np.sqrt(G*Mass/radii**2)*u.km/u.s #use formula to compute
        
        return VHern
   
    

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 14:43:17 2025

@author: lauramack
"""

# # In Class Lab 1
# Must be uploaded to your Github repository under a "Labs/Lab1" by midnight thursday

# ## Part A:  The Local Standard of Rest
# Proper motion of Sgr A* from Reid & Brunthaler 2004
# $\mu = 6.379$ mas/yr 
# 
# Peculiar motion of the sun, $v_\odot$ = 12.24 km/s  (Schonrich 2010)
# 
# 
# $v_{tan} = 4.74 \frac{\mu}{\rm mas/yr} \frac{R_o}{\rm kpc} = V_{LSR} + v_\odot$
# 
# 
# ### a)
# 
# Create a function called VLSR to compute the local standard of rest (V$_{LSR}$).
# 
# The function should take as input: the solar radius (R$_o$), the proper motion (mu)
# and the peculiar motion of the sun in the $v_\odot$ direction.
# 
# Compute V$_{LSR}$ using three different values R$_o$: 
# 1. Water Maser Distance for the Sun :  R$_o$ = 8.34 kpc   (Reid 2014 ApJ 783) 
# 2. GRAVITY Collaboration Distance for the Sun:  R$_o$ = 8.178 kpc   (Abuter+2019 A&A 625)
# 3. Value for Distance to Sun listed in Sparke & Gallagher : R$_o$ = 7.9 kpc 
# 



# Import Modules 
import numpy as np # import numpy
import astropy.units as u # import astropy units
from astropy import constants as const # import astropy constants

def VLSR(Ro, mu = 6.379, v_sun = 12.24*u.km/u.s):
    '''
    This function will comupte the velocity at the local standard of rest
    VLSR = 4.74*mu*Ro - vsun

    Parameters
    ----------
    Ro : astropy quantity
        distance from the sun to the galactic center (astropy units kpc)
    mu : astropy quanity, optional
        The proper motion of Sag A* (mas/yr). The default is 6.379 from 
        Reid & Brunthaler 2004
    v_sun : astropy quantity, optional
        the peculiar motion of the sun in the v direction (astropy units km/s). 
        The default is 12.24*u.km/u.s. from Schonrich+2010

    Returns
    -------
    VLSR : float
        the local standard of rest (astropy units km/s)

    '''
    VLSR = 4.74*mu*(Ro/u.kpc)*u.km/u.s - v_sun
    return VLSR

#differnt values of the distance to the galactic center
RoReid = 8.34*u.kpc # Reid + 2014
RoAbuter = 8.178*u.kpc #Abuter+2019
RoSparke = 7.9*u.kpc #Sparke and Gallagher text

#comput VSLR
VLSR_Reid = VLSR(RoReid)
print("Reid VLSR", np.round(VLSR_Reid, 3))

#VLSR from GRAVITY collaboration
VLSR_Abuter = VLSR(RoAbuter)
print("Abuter VLSR", np.round(VLSR_Abuter, 3))

#from Sparke and Gallagher
VLSR_Sparke = VLSR(RoSparke)
print("Sparke VLSR", np.round(VLSR_Sparke, 3))
print("\n")

# ### b)
# 
# compute the orbital period of the sun in Gyr using R$_o$ from the GRAVITY Collaboration (assume circular orbit)
# 
# Note that 1 km/s ~ 1kpc/Gyr

# orbitial period = 2piR/V
def TorbSun(Ro, Vc):
    '''
    A function that comuptes the orbital period of the sun
    T = 2 pi R / V

    Parameters
    ----------
    Ro : astropy quantity
        The distance to the galactic center from the sun (kpc)
    Vc : astropy quantity
        velocity of the sun the "v" direction (km/s)

    Returns
    -------
    T : astropy quantity
        Orbital period in Gyr

    '''
    VkpcGyr = Vc.to(u.kpc/u.Gyr) #conver V to kpc/Gyr
    T = 2*np.pi*Ro/VkpcGyr
    return T
    
VsunPec = 12.24*u.km/u.s #peculiar motion
Vsun = VLSR_Abuter + VsunPec # the total motion of the sun in v direction

##Orbital Period of sun
T_Abuter = TorbSun(RoAbuter, Vsun)
print("orbital period of the sun", np.round(T_Abuter, 3))
print("\n")

# ### c)
# 
# Compute the number of rotations about the GC over the age of the universe (13.8 Gyr)

AgeUniverse = 13.8*u.Gyr
print("number of times the sun has gone around the GC", np.round(AgeUniverse/T_Abuter, 5))
print("\n")


# ## Part B  Dark Matter Density Profiles
# 
# ### a)
# Try out Fitting Rotation Curves 
# [here](http://wittman.physics.ucdavis.edu/Animations/RotationCurve/GalacticRotation.html)
# 
# 
# ### b)
# 
# 
# In the Isothermal Sphere model, what is the mass enclosed within the solar radius (R$_o$) in units of M$_\odot$? 
# 
# Recall that for the Isothermal sphere :
# $\rho(r) = \frac{V_{LSR}^2}{4\pi G r^2}$
# 
# Where $G$ = 4.4985e-6 kpc$^3$/Gyr$^2$/M$_\odot$, r is in kpc and $V_{LSR}$ is in km/s
# 
# What about at 260 kpc (in units of  M$_\odot$) ? 

Grav = const.G.to(u.kpc**3/u.Gyr**2/u.Msun)

#Density Profile = VLSR^2/(4*pi*G*R^2)
#Mass(r) = integrate rho dV
#        = integrate rho 4*pi*r^2*dr
#        = integrate VLSR^2/G dr
#        = VLSR^2/G * r

def massIso(r, VLSR):
    '''
    This function will compute the dark matter mass enclosed withing a given 
    distance (r) assuming an iso thermal sphere model. M(r) = VLSR^2/G *r

    Parameters
    ----------
    r : astropy quantity
        distance from galactic center (kpc).
    VLSR : astropy quantity
        velocity at the local standard of rest (km/s)

    Returns
    -------
    M : astropy quantity
        mass enlosed withing r (Msun)

    '''
    VLSRkpcGyr = VLSR.to(u.kpc/u.Gyr) #translating to kpc/Gyr
    M = VLSRkpcGyr**2/Grav*r # iso thermal sphere mass profile
    return M

#mass enclosed with Ro (Gravity Colab)
mIsoSolar = massIso(RoAbuter, VLSR_Abuter)
print(f"{mIsoSolar:.2e}", "enclosed mass at solar radius")

mIso260 = massIso(260*u.kpc, VLSR_Abuter)
print(f"{mIso260:.2e}", "enclosed mass at 260 kpc")

# ## c) 
# 
# The Leo I satellite is one of the fastest moving satellite galaxies we know. 
# 
# 
# It is moving with 3D velocity of magnitude: Vtot = 196 km/s at a distance of 260 kpc (Sohn 2013 ApJ 768)
# 
# If we assume that Leo I is moving at the escape speed:
# 
# $v_{esc}^2 = 2|\Phi| = 2 \int G \frac{\rho(r)}{r}dV $ 
# 
# and assuming the Milky Way is well modeled by a Hernquist Sphere with a scale radius of $a$= 30 kpc, what is the minimum mass of the Milky Way (in units of M$_\odot$) ?  
# 
# How does this compare to estimates of the mass assuming the Isothermal Sphere model at 260 kpc (from your answer above)

# Potenial for Hernquist Sphere
#phi = -G*M/(r+a)
#escape speed becomes:
    #Vesc^2 2*G*M/(r+a)
#rearragne for M
# M = vesc^2/2/G*(r+a)

def massHernVesc(vesc, r, a=30*u.kpc):
    '''
    This function determines the total dar matter mass needed given an escape 
    speed, assuming a Hernquist profile
    
    Inputs:
        vesc(astropyquantity) escape speed (or speed of satellits) (km/s)
        r : (astropy quantity)  distance from the Galactic Center (kpc)
        a : (astropy quantity) the Herquist scale length (kpc)
        
    '''
    vescKpcGyr = vesc.to(u.kpc/u.Gyr) #translate to kpc/Gyr
    
    M = vescKpcGyr**2/2/Grav*(r+a)
    
    return M

Vleo = 196*u.km/u.s #speed of Leo 1
r = 260*u.kpc


MLeoI = massHernVesc(Vleo, r)
print(f"{MLeoI:.2e}", "mass at 260 kpc required to keep LeoI bound")



    




#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 28 10:09:38 2025

@author: lauramack
"""

import numpy as np
import astropy.units as u

def Read(filename):
    file = open(filename, 'r')
    
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    line2 = file.readline()
    label, value = line2.split()
    num_particles = float(value)
    
    file.close()
    
    data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    return time, num_particles, data
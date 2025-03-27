# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass2 import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit




class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """ initalize the class given a filename (ex: MW_000.txt) from the 
            simulation"""

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33COM = CenterOfMass("M33_000.txt", 2)
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33pos = M33COM.COM_P(0.1, 4)
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33vel = M33COM.COM_V(M33pos[0], M33pos[1], M33pos[2])
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass("M31_000.txt", 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31pos = M31COM.COM_P(0.1, 2)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31vel = M31COM.COM_V(M31pos[0], M31pos[1], M31pos[2])
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = (M33pos - M31pos).value
        self.v0 = (M33vel - M31vel).value
        
        ### get the mass of each component in M31 
        ### disk
        self.rdisk =  5 #scale length (kpc)
        
        #set with ComponentMass function. Remember to *1e12 to get the right 
        #units. Use the right ptype
        self.Mdisk  = ComponentMass("M31_000.txt", 2) * 1e12 #units of M_sun
        
        
        ### bulge
        self.rbulge = 1 #set scale length (kpc)

        self.Mbulge  = ComponentMass("M31_000.txt", 3) * 1e12 
        #set with ComponentMass function. Remember to *1e12 to get the right 
        #units Use the right ptype
        
        # Halo
        self.rhalo = 62  #set scale length from HW5 (kpc)

        self.Mhalo  = ComponentMass("M31_000.txt", 1) * 1e12
        #set with ComponentMass function. Remember to *1e12 to get the right 
        #units. Use the right ptype
     
    
    
    def HernquistAccel(self, M, r_a, r): # it is easiest if you take as an 
    #input the position VECTOR 
        """
        find the gravitational acceleration given by the halo or bulge 
        particles determined by a Hernquist profile
        
        Inputs: 
            M: Mass of halo or bulge (float)
            r_a: scale length of profile (float)
            r: position vector (array)
        Returns:
            Hern: acceleration from Hernquist profile (array)
        """

        ### **** Store the magnitude of the position vector
        rmag = (r[0]**2 + r[1]**2 + r[2]**2)**0.5
        
        ### *** Store the Acceleration
        Hern =  -(self.G * M)/(rmag * (r_a + rmag)**2) * r 
        #follow the formula in the HW instructions
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):
        # it is easiest if you take as an input a position VECTOR  r 
        '''
        Find the gravitation acceleration of the disk particles using a Miyamo 
        Nagai profile
        
        Inputs:
            M: disk mass (float)
            r_d: scale length (int) (should be self.rdisk)
            r: position vector for disk particles (array)
            
        Returns:
            a_MN: Miyamoto Magai acceleration (array)
        '''
        
        #a_MN = -GM/(R**2+B**2)**0.5 * (1, 1, B/(z**2+z_d**2)**0.5)
    
        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        R = np.sqrt(r[0]**2 + r[1]**2)
        
        z_d = self.rdisk/5
        B = r_d + np.sqrt(r[2]**2 + z_d**2)
        
        z_corr = B/((r[2]**2+z_d**2)**0.5) # the z direction is different from x,y
        mkvector = np.array([1,1,z_corr]) #make array to ensure that a vector is returned
        # notice z component is different
        
        
        a_MN = -(self.G*M)/ (R**2+B**2)**1.5 * r * mkvector
        
       
        return a_MN
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
            """ sums all acceleration components
            
            Inputs:
                r: position vector (array)
                
            Returns:
                accel: total aceleration for M31 (array)
            """

            ### Call the previous functions for the halo, bulge and disk
            # **** these functions will take as inputs variable we defined in the initialization of the class like 
            # self.rdisk etc. 
            Halo = self.HernquistAccel(self.Mhalo, self.rhalo, r)
            Bulge = self.HernquistAccel(self.Mbulge, self.rbulge, r)
            Disk = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)
            
            accel = Halo + Disk + Bulge
            
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
            return accel
    
    
    
    def LeapFrog(self, r, v, dt): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ update positions and velocties of M31 by guessing at next center 
        of mass position and updating from acceleration at the half timestep
        
        Inputs:
            r: position vector (array)
            v: velocity vector (array)
            dt: time interval (integer)
        
        Returns:
            rnew: the updated position vector (array)
            vnew: the updated velocity vector (array)"""
        
        # predict the position at the next half timestep
        rhalf = r + v * dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        ahalf = self.M31Accel(rhalf)
        vnew = v + ahalf*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew * dt/2
        
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ Find M31s orbit by looping over Leapfrog from start to end time
        
        Inputs: 
            t0: inital time (float)
            dt: timestep (float)
            tmax: maximum time (float)
            
        """

        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        rows = int(tmax/dt)+2
        orbit = np.zeros((rows, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        r = self.r0
        v = self.v0
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t<=tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            # **** store the new time in the first column of the ith row
            #orbit[i,0] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            
            r, v = self.LeapFrog(r, v, dt)
            orbit[i] = t, *tuple(r), *tuple(v)
            #r = orbit[i-1, 1:4]
            #v = orbit[i-1, 4:7]
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            #orbit[i, 1:4] = rnew
           
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
           # orbit[i, 4:7] = vnew
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
            
        #print(orbit[21])

        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#',
              #delimiter = " ", 
        header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
        .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

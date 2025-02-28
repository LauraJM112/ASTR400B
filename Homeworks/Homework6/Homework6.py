

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G


# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass2 import CenterOfMass




def OrbitCOM(galaxyname, start, end, n):
    """function that loops over all the desired snapshots to compute the COM 
        pos and vel as a function of time.
    inputs:
        galaxyname (string): the name of the galaxy in question (ex:MW)
        start (integer): the first snapshot to be calculated
        end (integer): the last snapshoot to be calculated
        n (integer): interval over which COM is calculated
          
    outputs: a file with the name "Orbit galaxyname.txt" with t,x,y,z,vx,vy,vz columns
    """
    
    
    # compose the filename for output
    fileout = "Orbit "+"%s"%(galaxyname)+".txt"
    print(fileout)

    
    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    delta = 0.1
    VolDec = 2
    if galaxyname == "M33": #set M33 value
        VolDec = 4
    
    
    # generate the snapshot id sequence 
    snap_id = np.arange(start, end, n)
    
    # it is always a good idea to also check if the input is eligible (not required)
    if len(snap_id) ==0:
        print("Invalid start, stop or n value") #should have code stop running here?
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((len(snap_id), 7))
    
    
    # a for loop 
    for  i in snap_id:# loop over files
        
        # compose the data filename (be careful about the folder)
        ilbl = '000' + str(i)
            # remove all but the last 3 digits
        ilbl = ilbl[-3:]
        filename="Homework6/"+"%s"%(galaxyname)+'/'+"%s"%(galaxyname) +"_"+ ilbl + '.txt'
        print(filename)
        # Initialize an instance of CenterOfMass class, using disk particles
        COM = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        PosCOM = COM.COM_P(delta, VolDec) #position
        x_COM = PosCOM[0] #grab x,y,z values
        y_COM = PosCOM[1]
        z_COM = PosCOM[2]
        
        VelCOM = COM.COM_V(x_COM, y_COM, z_COM) #velocity
        vx_COM = VelCOM[0] #grab vx,vy,vz values
        vy_COM = VelCOM[1]
        vz_COM = VelCOM[2]
    
        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        time = COM.time/1000 #make gigayears
        values = [time.value, x_COM.value, y_COM.value, z_COM.value, 
                  vx_COM.value, vy_COM.value, vz_COM.value] #list of values to put into orbit array
        index = np.where(snap_id == i) #get index for snap_id and orbit
        orbit[index] = values #put values into orbit
        

        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    
    np.savetxt(fileout, orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        




# Recover the orbits and generate the COM files for each galaxy
# read in 800 snapshots in intervals of n=5
# Note: This might take a little while - test your code with a smaller number of snapshots first! 
#OrbitCOM('MW', 0, 800, 5)
#OrbitCOM("M33", 0, 800, 5)
#OrbitCOM("M31", 0, 800, 5)



# Read in the data files for the orbits of each galaxy that you just created
# headers:  t, x, y, z, vx, vy, vz
# using np.genfromtxt
#Milky Way
dataMW = np.genfromtxt("Orbit MW.txt")
MWtime = dataMW[:,0]
MWx = dataMW[:,1]
MWy = dataMW[:,2]
MWz = dataMW[:,3]
MWvx = dataMW[:,4]
MWvy = dataMW[:,5]
MWvz = dataMW[:,6]
MWPosition = [MWx, MWy, MWz]
MWVelocity = [MWvx, MWvy, MWvz]


#M33
dataM33 = np.genfromtxt("Orbit M33.txt")
M33time = dataM33[:,0]
M33x = dataM33[:,1]
M33y = dataM33[:,2]
M33z = dataM33[:,3]
M33vx = dataM33[:,4]
M33vy = dataM33[:,5]
M33vz = dataM33[:,6]
M33Position = [M33x, M33y, M33z]
M33Velocity = [M33vx, M33vy, M33vz]
#print(M33x)

#M31
dataM31 = np.genfromtxt("Orbit M31.txt")
M31time = dataM31[:,0]
M31x = dataM31[:,1]
M31y = dataM31[:,2]
M31z = dataM31[:,3]
M31vx = dataM31[:,4]
M31vy = dataM31[:,5]
M31vz = dataM31[:,6]
M31Position = np.array([M31x, M31y, M31z])
M31Velocity = np.array([M31vx, M31vy, M31vz])



# function to compute the magnitude of the difference between two vectors 
# You can use this function to return both the relative position and relative velocity for two 
# galaxies over the entire orbit  

def VectorDiff (V1, V2):
    """
    Caluculates the diffrnce between to 3D vectors and its magnitude

    Parameters
    ----------
    V1 : 3D numpy array of floats
        Vector 1
    V2 : 3D numpy array of floats
        Vector 3

    Returns
    -------
    mag : array of floats
        contains the magnitude of diffrence between the vectors
    """

    Diff = np.subtract(V1,V2)
    mag = (Diff[0]**2+Diff[1]**2+Diff[2]**2)**0.5
    
    return mag
    


# Determine the magnitude of the relative position and velocities 

# of MW and M31
MW_M31Pos = VectorDiff(MWPosition, M31Position)
#print(MW_M31Pos)

MW_M31Vel = VectorDiff(MWVelocity, M31Velocity)



# of M33 and M31
M31_M33Pos = VectorDiff(M31Position, M33Position)
#print(M31_M33)
M31_M33Vel = VectorDiff(M31Velocity, M33Velocity)


# Plot the Orbit of the galaxies 
#################################

plt.plot(MWtime, MW_M31Pos)
plt.title("MW and M31 Separation as a function of time")
plt.xlabel("Time (Gyr)")
plt.ylabel("Separation (kpc)")
plt.savefig("MW-M31 Separation")
plt.show()

plt.plot(MWtime, M31_M33Pos)
plt.title("M31 and M33 Separation as a function of time")
plt.xlabel("Time (Gyr)")
plt.ylabel("Separation (kpc)")
plt.savefig("M31-M33 Separation")
plt.show()


# Plot the orbital velocities of the galaxies 
#################################
plt.plot(MWtime, MW_M31Vel)
plt.title("MW and M31 relative velocities as a function of time")
plt.xlabel("Time (Gyr)")
plt.ylabel("Velocity (kpc/s)")
plt.savefig("MW-M31 Relative Velocity")
plt.show()

plt.plot(MWtime, M31_M33Vel)
plt.title("M31 and M33 relative velocities as a function of time")
plt.xlabel("Time (Gyr)")
plt.ylabel("Velocity (kpc/s)")
plt.savefig("M31-M33 Relative Velocity")
plt.show()


# Overplot Seperation and Relative Velocity
# MW and M#1
plt.plot(MWtime, MW_M31Pos, label = "position")
plt.plot(MWtime, MW_M31Vel, label = "relative velocity")
plt.title("MW and M31")
plt.legend()
plt.show()

#M31 and M33
plt.plot(MWtime, M31_M33Pos, label = "Position")
plt.plot(MWtime, M31_M33Vel, label = "relative velocity")
plt.title("M31 and M33")
plt.legend()
plt.show()

#What happens when merger happens
#MW and M31 Merge (zoom in on plot)
plt.plot(MWtime, MW_M31Pos)
plt.title("MW and M31 Separation as a function of time")
plt.xlabel("Time (Gyr)")
plt.ylabel("Separation (kpc)")
#plt.yscale("log")
plt.xlim(5.75)
plt.show()


#what happens to M33 (zoom in)
plt.plot(MWtime, M31_M33Pos, label = "Position")
plt.title("M33 an M31 Positions")
plt.axvline(x=6.75, c = 'r', label = "time of merger", linewidth=0.75)
#plt.yscale('log')
plt.legend()
#plt.xlim(7.5)
plt.show()


#Q5 - 
plt.plot(MWtime, M31_M33Pos, label = "Position")
plt.title("M33 an M31 Positions")
plt.axvline(x=6.75, c = 'r', label = "time of merger", linewidth=0.75)
#plt.yscale('log')
plt.scatter(7.5, 109, label = "apocenter (t=7.5 , sep = 109)")
plt.scatter(8.9, 90, label = "aopcenter (t=8.9, sep=90)")
plt.legend()
#plt.xlim(7.5)
plt.show()


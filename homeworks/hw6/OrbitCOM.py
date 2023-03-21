

# Homework 6 Template
# G. Besla & R. Li




# import modules
import numpy as np
import astropy.units as u
from astropy.constants import G

# import plotting modules
import matplotlib.pyplot as plt

# my modules
from readfile import Read
# Step 1: modify CenterOfMass so that COM_P now takes a parameter specifying 
# by how much to decrease RMAX instead of a factor of 2
from CenterOfMass import CenterOfMass



def OrbitCOM(galaxy, start, end, n):
    """function that loops over all the desired snapshots to compute the COM pos and vel as a function of time.

    inputs:
        galaxy: string
            Name of the galaxy, e.g. 'MW'
        start: int
            Number of the first snapshot
         end: int
            Number of the last snapshot
        n: int
            Number of intervals over which the COM will be returned. 

    outputs: 
    """
    
    # compose the filename for output
    
    outfile = f"Orbit-{galaxy}.txt"

    #  set tolerance and VolDec for calculating COM_P in CenterOfMass
    # for M33 that is stripped more, use different values for VolDec
    
    # generate the snapshot id sequence 
    # it is always a good idea to also check if the input is eligible (not required)
    assert type(start) == int
    assert type(end) == int
    assert start - end != 0
    
    snap_ids = np.arange(start, end, n) 
    
    # initialize the array for orbital info: t, x, y, z, vx, vy, vz of COM
    orbit = np.zeros((len(snap_ids), 7))
    
    # a for loop 
    for i, snap_id in enumerate(snap_ids): # loop over files
        
        # compose the data filename (be careful about the folder)
        filename = f'data/{galaxy}_{["000"+str(snap_id)][0][-3:]}.txt'
        # Initialize an instance of CenterOfMass class, using disk particles
        M = CenterOfMass(filename, 2)
        # Store the COM pos and vel. Remember that now COM_P required VolDec
        if galaxy == "M33":
            M_COM_P = M.COM_P(volDec=4)
        else:
            M_COM_P = M.COM_P()

        M_COM_V = M.COM_V(M_COM_P[0],M_COM_P[1],M_COM_P[2])

        # store the time, pos, vel in ith element of the orbit array,  without units (.value) 
        # note that you can store 
        # a[i] = var1, *tuple(array1)
        orbit[i] = M.time.value / 1000, *tuple(M_COM_P.value), *tuple(M_COM_V.value)
        
        # print snap_id to see the progress
        print(snap_id)
        
    # write the data to a file
    # we do this because we don't want to have to repeat this process 
    # this code should only have to be called once per galaxy.
    np.savetxt(f'output/{outfile}', orbit, fmt = "%11.3f"*7, comments='#',
               header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                      .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))

def vector_diff(a, b):
    """ Computes magnitude of the difference between two vectors.
    
        Inputs:
            a: numpy array
                First Vector (x, y, z)
            b: numpy array
                Second Vector (x, y, z)

        Outputs:
            mag: float
                Magnitude of the difference vector.
    """
    sub = a - b
    res = np.sqrt(sub[0]**2 + sub[1]**2 + sub[2]**2)
    return np.abs(res)


def main():
    ## as given in the instructions, we will be going from snapshot 0 to 800
    ## with n=5. I will set this up to run on nimpy and then just leave it
    ## running.
    ## for i in ['MW', 'M31', 'M33']:
    ##     OrbitCOM(i, 0, 800, 5)

    # Recover the orbits and generate the COM files for each galaxy
    # read in 800 snapshots in intervals of n=5
    # Note: This might take a little while - test your code with a smaller number of snapshots first! 
    
    
    
    
    # Read in the data files for the orbits of each galaxy that you just created
    # headers:  t, x, y, z, vx, vy, vz
    # using np.genfromtxt
    
    MW = np.genfromtxt("output/LowRes/Orbit-MW.txt")
    M31 = np.genfromtxt("output/LowRes/Orbit-M31.txt")
    M33 = np.genfromtxt("output/LowRes/Orbit-M33.txt")

    # function to compute the magnitude of the difference between two vectors 
    # You can use this function to return both the relative position and relative velocity for two 
    # galaxies over the entire orbit  
    
    time = MW[:, 0]

    ## MW vs M31
    
    # Determine the magnitude of the relative position and velocities 
    
    # of MW and M31
    MW_M31_Separation = vector_diff(np.array([MW[:, 1], MW[:, 2], MW[:, 3]]), np.array([M31[:, 1], M31[:, 2], M31[:, 3]],))
    MW_M31_Velocity = vector_diff(np.array([MW[:, 4], MW[:, 5], MW[:, 6]]), np.array([M31[:, 4], M31[:, 5], M31[:, 6]],))

    # of M33 and M31
    M33_M31_Separation = vector_diff(np.array([M33[:, 1], M33[:, 2], M33[:, 3]]), np.array([M31[:, 1], M31[:, 2], M31[:, 3]],))
    M33_M31_Velocity = vector_diff(np.array([M33[:, 4], M33[:, 5], M33[:, 6]]), np.array([M31[:, 4], M31[:, 5], M31[:, 6]],))

    out = np.array([time, M33_M31_Separation, M33_M31_Velocity])
    print(out)
    np.savetxt("M31_M33_orbit.txt", out, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
    fig, ax = plt.subplots(2, 2, sharey="row", sharex=True, figsize=(8, 4))
    
    # Plot the Orbit of the galaxies 
    #################################
    ax[0, 0].plot(time, MW_M31_Separation)
    ax[0, 0].set_title("MW vs M31 (Separation)")

    ax[0, 1].plot(time, M33_M31_Separation)
    ax[0, 1].set_title("M33 vs M31 (Separation)")

    # Plot the orbital velocities of the galaxies
    #################################
    ax[1, 0].plot(time, MW_M31_Velocity)
    ax[1, 0].set_title("MW vs M31 (Velocity)")

    ax[1, 1].plot(time, M33_M31_Velocity)
    ax[1, 1].set_title("M33 vs M31 (Velocity)")

    ax[1, 0].set_xlabel("Time [Gyr]")
    ax[0, 0].set_ylabel("Separation [kpc]")
    ax[1, 0].set_ylabel("V [km/s]")

    plt.savefig("fig.png")


if __name__ == "__main__":
    main()

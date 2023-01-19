import numpy as np
import astropy.units as u
import sys

def Read(name, ):
    """ 
    Reads data from a tab-separated text file called "MW_000.txt containing galactic simulation data. Returns the time, total no. of particles, and all the data within.
    Inputs : filename
    Returns : time in timestep units, no. of particles in the file, and all remaining filedata in a numpy array.
    """
    with open(name, 'r') as file:
        time = float(file.readline().split()[1]) * u.Myr ## time looks like a float in the file. should it be an int instead?
        particle_no = file.readline().split()[1]
    data = np.genfromtxt(name, dtype=None, names=True, skip_header=3) ## extract data from file
    return time, particle_no, data

def main():
    Read("MW_000.txt", )

if __name__ == "__main__":
    main()

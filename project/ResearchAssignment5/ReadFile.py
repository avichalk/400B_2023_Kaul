# imports
import numpy as np
import astropy.units as u

# ---------- FUNCTIONS ----------- #

def Read(filename):
	'''
	Reads ASTR400B simulation snapshots. 

	PARAMETERS
	---------
	filename : `~str`
		Name of the snapshot to be read

	RETURNS
	-------
	time : `astropy.quantity`
		Time of snapshot in Myr
	N : `int`
		Number of particles in the snapshot
	data : `numpy.ndarray`
		Particle data, contains particle type, mass (in 1e10 Msun), x, y, z positions in kpc, and vx, vy, vz 
		velocities in km/s
	'''

	file = open(filename, 'r') # open file

	# read and store time 
	line1 = file.readline()
	label, value = line1.split()
	time = float(value)*u.Myr

	# read and store particle number
	line2 = file.readline()
	label, value = line2.split()
	N = int(value)

	# read and store data frame with particle data
	data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

	# return outputs
	return time, N, data

# --------- MAIN ---------- # 
if __name__ == "__main__":

	time, N, data = Read("MW_020.txt") # read snap
	# print checks
	print("time =", time)
	print("N =", N)
	print("Mass of first particle =", data['m'][0]*u.Msun*1e10)
	print("Type of first particle =", data['type'][0])
	

# imports
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from ReadFile import Read

# ------------- FUNCTIONS -------------- #
def plot(filename, partType):
	'''
	Returns magnitude of the position and velocity vectors and the mass for specified particles.

	PARAMETERS
	----------
	filename : `str`
		Name of the snapshot to be read
	partType : `int`
		Type of particle
	partNum : `int or array-like`
		Particle (line) numbers we want data for

	RETURNS
	-------
	r : `astropy.quantity`
		3-D distance from simulation origin, in kpc
	v : `astropy.quantity`
		Mag of 3-D velocity in km/s
	m : `astropy.quantity`
		Particle mass in Msun
	'''

	# read particle data from snap
	t, N, data = Read(filename)

	# select only particles of desired type
	idx = np.where(data['type'] == partType)

	# make arrays of particle data for the desired type and line number
	m = data['m']#[idx]*1e10*u.Msun # also convert to Msun
	x = data['x']#[idx]*u.kpc
	y = data['y']#[idx]*u.kpc
	z = data['z']#[idx]*u.kpc
	vx = data['vx']#[idx]*u.km/u.s
	vy = data['vy']#[idx]*u.km/u.s
	vz = data['vz']#[idx]*u.km/u.s
	print(x, y, z)

	fig = plt.figure()
	ax = fig.add_subplot(projection='3d')
	ax.scatter(x, y, z)
	plt.show()

# ---------- MAIN ----------- #
if __name__ == "__main__":

	plot("MW_000.txt", 1)
	
	
	

## The goal of my research project is to track the mass transfer between M31 and the Milky Way using Simulation Data.
## For this assignment, we will only concern ourselves with the position of the transferred mass, as opposed to the kinematics.
## We will be creating a plot showing the mass transferred between the two galaxies at different snapshots.


import numpy as np # type: ignore
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt

from readfile import Read
from CenterOfMass import CenterOfMass
from mass_distribution import MassProfile
from astropy.constants import G
G = G.to(u.kpc * u.km**2 / u.s**2 / u.M_sun)

class GalaxyPos:
    def __init__(self, galaxy_name):
        self.time, self.total, self.data = Read(galaxy_name)                                                                                             
        print(self.data.dtype)

    def position(self, ptype):
        index = np.where(self.data['type'] == ptype)
        
        self.x = self.data['x'][index]
        self.y = self.data['y'][index]
        self.z = self.data['z'][index]

        return self.x, self.y

    def potential(self, ptype):
        ## i understand fitting to a Hernquist profile, but how to find the scale radius?
        pass
"""
At each timestep, calculate all particle potential energy and kinetic energy. use hernquist profile and COM to find which particles are bound to what. 
So, at each timestep, a new class is instantiated.
"""

class GenerateEnergies:
    def __init__(self, snapshot, galaxy):
        """
        Generate energies for each particle and save back into a file
        that can be read from for future work.

        We calculate the binding for all particles of one galaxy first, then all particles of the second galaxy.
        """
        self.name = self.filename(snapshot, galaxy) 
        self.shot = snapshot
        self.galaxy = galaxy

        self.time, self.total, self.data = Read(self.name)

        self.m = np.transpose(self.data['m'])
        self.x = np.transpose(self.data['x']) #* u.kpc
        self.y = np.transpose(self.data['y']) #* u.kpc
        self.z = np.transpose(self.data['z']) #* u.kpc
        
        self.vx = np.transpose(self.data['vx'])
        self.vy = np.transpose(self.data['vy'])
        self.vz = np.transpose(self.data['vz'])


    def filename(self, snapshot, galaxy):
        return f'{galaxy}_{["000"+str(snapshot)][0][-3:]}.txt'

    def KE(self,):
        """
        Generates KE for each particle in the dataset.
        KE = 1/2mv^2. Returns a sequential array.
        Inputs:
            self.data : `np array of floats`
                Contains m : mass in 1e10 solar masses
                vx, vy, vz : component velocities in km/s
        """
        print(f"Calculating KE.")
        v_squared = self.vx**2 + self.vy**2 + self.vz**2 
        return 0.5 * self.m * v_squared * (u.km / u.s)**2 ## * 1e10 * u.Msun * 

    def V(self,):
        """
        Generates V for each particle in the dataset.
        We find the mass enclosed using a Hernquist profile.
        Unclear ATM how accurate this is.
        Inputs:
            self.data : `np array of floats`
                Contains m : mass in 1e10 solar masses
                vx, vy, vz : component velocities in km/s
        Calculates separately for MW and M31, and returns in array.
        """
        ## permutations of measuring V?
        ## this entire function to be called twice for both data (MW and M31)
        out = []
        for i in {"MW", "M31"}:
            print(f"Calculating V for {i}.")
            COM = CenterOfMass(self.filename(self.shot, i), 2)
            COM_p = COM.COM_P(0.1).value
            a = self.x - COM_p[0]
            b = self.y - COM_p[1]
            c = self.z - COM_p[2]
            r = np.sqrt(np.square(a)+np.square(b)+np.square(c))
            Profile = MassProfile(i, self.shot) ## ansc inputs, don't actually work
            print("Calculating mass enclosed")
            M = Profile.MassEnclosedTotal(r)
            print("doing V")
            V = - G * M / r
            np.savetxt(f"{i}", V)
            out.append(V)
        return out

    def check_if_bound(self,):
        """
        Check total energy. If -ve, it is gravitationally bound to that galaxy.
        If positive, unbound. Check against both galaxies.
        """
        KE = self.KE()
        MW_V, M31_V = self.V()
        
        print(f"Calculating Total E.")
        E_MW = KE.value + MW_V.value
        E_M31 = KE.value + M31_V.value

        return np.sign(E_MW), np.sign(E_M31)



def main():
    ## i'm going to plot each different particle type in its own plot
    #MW = GalaxyPos("MW_000.txt")
    #M31 = GalaxyPos("M31_000.txt")
    #
    ### Plot COM so we can see how it changes over time. We can also use this to fit a Hernquist profil
    #MW_COM = CenterOfMass("MW_000.txt", 2)
    ##MW_COM_p = MW_COM.COM_P()

    ### loop over all particles and check a) unbound to host galaxy b) bound to new galaxy. if neither are true, assume still bound to host galaxy
    #
    #fig, ax = plt.subplots(figsize=[10, 5])

    #ax.scatter(*MW.position(2),label="MW")
    #ax.scatter(*M31.position(2),label="M31")

    #ax.legend()
    #plt.show()

    for i in {"MW", "M31"}:
        E = GenerateEnergies(0, i)
        res = E.check_if_bound()
        print(res)


if __name__ == "__main__":
    main()

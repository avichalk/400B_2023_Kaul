## The goal of my research project is to track the mass transfer between M31 and the Milky Way using Simulation Data.
## For this assignment, we will only concern ourselves with the position of the transferred mass, as opposed to the kinematics.
## We will be creating a plot showing the mass transferred between the two galaxies at different snapshots.

import numpy as np
import astropy.constants as const
import matplotlib.pyplot as plt

from readfile import Read

class GalaxyPos:
    def __init__(self, galaxy_name):
        self.time, self.total, self.data = Read(galaxy_name)                                                                                             

    def position(self, ptype):

        index = np.where(self.data['type'] == ptype)
        
        self.x = self.data['x'][index]
        self.y = self.data['y'][index]
        self.z = self.data['z'][index]

        return self.x, self.y

def main():
    ## i'm going to plot each different particle type in its own plot
    filename = lambda galaxy, snap_id: f'../hw6/data/LowRes/{galaxy}/{galaxy}_{["000"+str(snap_id)][0][-3:]}.txt'

    for i in range(0, 500, 50):
        MW = GalaxyPos(filename("MW", i))
        M31 = GalaxyPos(filename("M31", i))
        
        fig, ax = plt.subplots(figsize=(8, 4))

        ax.scatter(*MW.position(2),label="MW")
        ax.scatter(*M31.position(2),label="M31")

        ax.legend()

        plt.savefig(f"output/{i}.png")
        

        print(i)


    ## issue : how to classify whether or not one part of the galaxy is inside another part? Use COM function?
    ## which center of mass is the galaxy closer to
    ## or something along the lines of an energy argument
    ## use half-mass radius, 
    ## check energy. calculate potential by getting potential of the particle against the galaxy (assume point mass), and you know the KE.
    ## find the total E, if it is -ve, it is gravitationally bound to that galaxy

    ## use binning algorithm to clear up graph. contour plots within galaxy?


if __name__ == "__main__":
    main()

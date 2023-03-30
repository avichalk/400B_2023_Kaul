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
    MW = GalaxyPos("MW_000.txt")
    M31 = GalaxyPos("M31_000.txt")
    
    fig, ax = plt.subplots()

    ax.scatter(*MW.position(2),label="MW")
    ax.scatter(*M31.position(2),label="M31")

    ax.legend()
    plt.show()

if __name__ == "__main__":
    main()

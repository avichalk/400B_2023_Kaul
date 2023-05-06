## Galactic Mergers are a process during which two or more galaxies fall together
## to form one galaxy. The structure of the galaxy undergoes massive changes, which 
## can have huge consequences for stellar formation and galactic evolution. Our aim 
## is to study the nature of the mass transfer during the merger between M31 and the
## MW. This will help further our understanding of galactic evolution, and help 
## explain where mass transferred between the galaxies ends up.


## This code has multiple classes. Towards the top, GalaxyPos is useful for plotting the
## position and velocity of the particles in the galaxy. Below, GenerateEnergies is used
## to generate the total energy of a particle within the system and save it. 
## And finally, we have the main function, which contains all the plotting code.

import numpy as np # type: ignore
import astropy.constants as const
import astropy.units as u
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys
import os
import gc ## garbage collection
from scipy.ndimage import gaussian_filter1d

from matplotlib.gridspec import GridSpec

import matplotlib.colors as colors
from matplotlib.colors import LogNorm

from ReadFile import Read
from CenterOfMass import CenterOfMass
from mass_distribution import MassProfile
from astropy.constants import G
G = G.to(u.kpc * u.km**2 / u.s**2 / u.M_sun)

def filename(snapshot, galaxy):
    """
    Returns the snapshot name on my local computer

    Inputs:
        snapshot : `int`
            Snapshot Number.
        Galaxy : `string`
            Galaxy name.

    Outputs:
        `string`
            The path to the snapshot and its name.
    """
    return f'snapshots/{galaxy}/{galaxy}_{["000"+str(snapshot)][0][-3:]}.txt'

class GalaxyPos:
    def __init__(self, snapshot, galaxy):
        """"
        Inputs:
            snapshot : `int`
                Snapshot Number.
            Galaxy : `string`
                Galaxy name
        """
        self.shot = snapshot
        self.galaxy = galaxy
        self.time, self.total, self.data = Read(filename(snapshot, galaxy))                                                                                             
        ## reads in all our data

    def plot_position(self, ptype=None):
        """
        Plots position of desired particles.
        Inputs:
            ptype : `int`
                Type of particle you want plotted
        """
        if ptype == None:
            self.x = self.data['x']
            self.y = self.data['y']
            self.z = self.data['z']
        else:
            index = np.where(self.data['type'] == ptype)
            
            self.x = self.data['x'][index]
            self.y = self.data['y'][index]
            self.z = self.data['z'][index]
        #print(self.x, self.y, self.z)
        
        fig = plt.figure()#, ax = plt.subplots(figsize=[10, 5])
        ax = fig.add_subplot(projection='3d')

        #col1 = np.where(mw_to_mw < mw_to_m31, 'r', np.where(mw_to_m31 < mw_to_mw, 'b', ''))
        #col2 = np.where(m31_to_mw < 0, 'r', np.where(m31_to_m31 < 0, 'b', ''))
        #m31pos = M31.position(2)    

        ax.scatter(self.x, self.y, self.z,  label="MW",) #marker='o') #c=col1,
        #ax.scatter(*m31pos, c=col2, label="M31")

        ax.legend()
        #plt.savefig(f"output/{i}.png")
        plt.show()

    def position(self, ptype=None, COM=False):
        """
        Returns position of all particles in Cartesian coordinates
        cooresponding to a particle particle type. COM logic was added
        for functionality where you could calculate the distance from
        the other galaxy, but I ran out of time :(

        Inputs:
            ptype : `None, int`
                Particle type. If none, we sum over all particle types.
            COM : `bool`
                If True, subtract COM. If false, do not.
        """
        if ptype == None:
            x = self.data['x']
            y = self.data['y']
            z = self.data['z']
        else:
            index = np.where(self.data['type'] == ptype)
            
            x = self.data['x'][index]
            y = self.data['y'][index]
            z = self.data['z'][index]

        if COM:
            range = [1, 2, 3] if ptype == None else [ptype]
            C = np.sum([CenterOfMass(filename(self.shot, self.galaxy), i) for i in range])
            COM_p = C.COM_P(0.1).value
            return x - COM_p[0], y - COM_p[1], z - COM_p[2]
        return x, y, z

    def velocity(self, ptype=None, COM=False):
        """
        Returns velocity components of all particles in Cartesian coordinates
        cooresponding to a particle particle type. COM logic was added
        for functionality where you could subtract COM velocity of
        the other galaxy, but I ran out of time :(

        Inputs:
            ptype : `None, int`
                Particle type. If none, we sum over all particle types.
            COM : `bool`
                If True, subtract COM. If false, do not.
        """

        if ptype == None:
            x = self.data['vx']
            y = self.data['vy']
            z = self.data['vz']
        else:
            index = np.where(self.data['type'] == ptype)
            
            x = self.data['vx'][index]
            y = self.data['vy'][index]
            z = self.data['vz'][index]

        if COM:
            range = [1, 2, 3] if ptype == None else [ptype]
            C = np.sum([CenterOfMass(filename(self.shot, self.galaxy), i) for i in range])
            COM_p = C.COM_P(0.1).value
            COM_v = C.COM_V(*COM_p).value
            return x - COM_v[0], y - COM_v[1], z - COM_v[2]

        return x, y, z

class GenerateEnergies:
    def __init__(self, snapshot, galaxy, ptype=None):
        """
        Generate energies for each particle and save back into a file
        that can be read from for future work.

        We calculate the binding for all particles of one galaxy first,
        then all particles of the second galaxy. The goal is to check if
        one is greater than the other.

        Inputs:
            snapshot : `int`
                Snapshot number
            galaxy : `string`
                Galaxy name
            ptype : `None, int`
                Particle type. If none, we sum over all particle types.
        """
        print(f"    Galaxy: {galaxy}")
        print(f"    Particle Type: {ptype}")
        self.name = filename(snapshot, galaxy) 
        self.shot = snapshot
        self.galaxy = galaxy
        self.ptype = ptype

        if len(sys.argv) > 1 and sys.argv[1] == "-r":
            ## flag to generate data. if not, we load from file.

            self.time, self.total, self.data = Read(self.name)

            
            if ptype == None:

                self.m = np.transpose(self.data['m'])
                self.x = np.transpose(self.data['x']) #* u.kpc
                self.y = np.transpose(self.data['y']) #* u.kpc
                self.z = np.transpose(self.data['z']) #* u.kpc
                
                self.vx = np.transpose(self.data['vx'])
                self.vy = np.transpose(self.data['vy'])
                self.vz = np.transpose(self.data['vz'])
            else:

                index = np.where(self.data['type'] == ptype)

                self.m = np.transpose(self.data['m'][index])
                self.x = np.transpose(self.data['x'][index]) #* u.kpc
                self.y = np.transpose(self.data['y'][index]) #* u.kpc
                self.z = np.transpose(self.data['z'][index]) #* u.kpc
                
                self.vx = np.transpose(self.data['vx'][index])
                self.vy = np.transpose(self.data['vy'][index])
                self.vz = np.transpose(self.data['vz'][index])
        else:

            print(f"    Forego loading data")



    def KE(self, COM):
        """
        Generates KE for each particle in the dataset.
        KE = 1/2mv^2. Subtracts COM velocity and returns a sequential array.

        Inputs:
            self.data : `np array of floats`
                Contains m : mass in 1e10 solar masses
                x, y, z : radii in kpc
                vx, vy, vz : component velocities in km/s
            COM : `class`
                COM object to calculate COM velocity

        Outputs:
            KE: `3d np array of floats`
                Kinetic energy in units of (km/s)^2 * 1e10 Msun
        """
        print(f"        Calculating KE.")
        
        vc = COM.COM_V(self.x * u.kpc, self.y * u.kpc, self.z * u.kpc) ## returns values in km/s
        vx, vy, vz = self.vx - vc[0].value, self.vy - vc[1].value, self.vz - vc[2].value
        v_squared = vx**2 + vy**2 + vz**2 ## subtract COM velocity
        return 0.5 * self.m * v_squared #* (u.km / u.s)**2 ## * 1e10 * u.Msun * 

    def V(self,):
        """
        Generates V for each particle in the dataset.
        We find the mass enclosed using a Hernquist profile.
        This class is run for one galaxy at a time (parent galaxy). So, it first calculates
        the potentials for all MW partiles against MW and M31, then for all
        M31 particles against MW and M31. 
        Inputs:
            self.data : `np array of floats`
                Contains m : mass in 1e10 solar masses
                x, y, z : radii in kpc
                vx, vy, vz : component velocities in km/s
        Outputs:
            out : `array of numpy arrays`
                Potential of parent galaxy's particles to (a) themselves and 
                (b) the other galaxy.
        """
        out = []

        print(f"        Calculating V for {self.galaxy}.") 
        
        ## now calculates potential for each particle for MW, then M31
        for i in {"MW", "M31"}:
            range = [1, 2, 3] if self.ptype == None else [self.ptype] 
            ## calculates distance from parent to COM
            COM = np.sum([CenterOfMass(filename(self.shot, i), j) for j in range])
            COM_p = COM.COM_P(0.1).value
            a = self.x - COM_p[0]
            b = self.y - COM_p[1]
            c = self.z - COM_p[2]
            r = np.sqrt(np.square(a)+np.square(b)+np.square(c))

            ## calculates mass profile
            print(f"        Calculating against {i}")
            Profile = MassProfile(i, self.shot)

            ## finds mass enclosed
            print("             Calculating mass enclosed")
            full_halo_mass_index = np.where(Profile.data["type"] == 1)[0]
            full_halo_mass = np.sum(Profile.data["m"][full_halo_mass_index]) ## scalar

            ## plugs into our equation to get the potential
            M = Profile.HernquistMass(r, 62, full_halo_mass)
            print("             Calculating V")
            V = - G * M / r
            out.append(V.value)

        return out

    def check_if_bound(self,):
        """
        Generates total energy, with some logic about when to generate
        and when to load from txt files instead. Runs very slowly...

        Inputs:
            N/A
        Outputs:
            E_MW, E_M31 : `numpy arrays of floats`
                Energy of parent galaxy against MW,
                Energy of parent galaxy against M31
        """
        
        if len(sys.argv) > 1 and sys.argv[1] == "-r":

            range = [1, 2, 3] if self.ptype == None else [self.ptype] 
            COM = np.sum([CenterOfMass(filename(self.shot, self.galaxy), i) for i in range])

            KE = self.KE(COM) ## units of (km/s)^2 * 1e10 Msun

            res = self.V() ## units of (km/s)^2 * 1e10 Msun
            np.savetxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}0.txt", res[0])
            np.savetxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}1.txt", res[1])
            np.savetxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}2.txt", KE)
            

        else:

            res = [0, 0]
            res[0] = np.loadtxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}0.txt")
            res[1] = np.loadtxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}1.txt")
            KE = np.loadtxt(f"potential/{self.galaxy}_{self.shot}_{self.ptype}2.txt")

        print(f"        Calculating Total E.")
        E_MW = KE + res[0]
        E_M31 = KE + res[1] ## unit errors?

        return (E_MW), (E_M31)



def main():
    """
    NOTE: Unfortunately, this function is a mess. The procedure I originally followed is that
    I would write some code, generate the requisite plots, and then comment the code out to write new code.
    Unfortunately, this has... caused some issues while going back and reconstructing what the code originally did.

    I'm almost certain that it works fine, now. The clunky way I've gotten around not generating new plots each time
    is to add a do_flag that can be changed to True if you want that function run. This is likely not best practise, and yes,
    I probably should have separated these into their own functions. Live and learn.

    General notes: The try, except blocks are there because some snapshots would randomly cause the code to error out.
    """
    ptypelist = [1, 2, 3]
    ptypecheatsheet = {1:"Halo", 2:"Disk", 3:"Bulge"}


    do_flag = False
    if do_flag:
         for i in range(300, 500, 20):
             print(f"Snapshot {i}.")
             fig, ax = plt.subplots(2, 3, figsize=[30, 15])
             for index, galaxy_name in enumerate(['M31', 'MW',], 0): 

                print(f"    {galaxy_name}")
                for idx, ptype in enumerate(ptypelist):
                    M = GalaxyPos(i, galaxy_name)
                    vx, vy, vz = M.velocity(ptype, COM=True)
                    vtot = np.sqrt(vx**2 + vy**2 + vz**2)

                    nbins = 500
                    ax[index, idx].hist(vtot, bins=nbins)
                    ax[index, idx].set_xlabel("velocity [km/s]")
                    ax[index, idx].set_ylabel("No. of stars at velocity")
                    ax[index, idx].set_title(f"Velocity distribution of {ptypecheatsheet[ptype]} particles in {galaxy_name}")
                fig.savefig(f"output/vhist/vhist{i}")
                plt.close()


    do_flag = False
    if do_flag:
        fig, ax = plt.subplots(2, ncols=len(ptypelist), figsize=[30, 15],) #squeeze=False) #subplot_kw=dict(projection='3d'))

        ## this code finds the mean and stdev of the galactocentric radius

        # I was coming up on the deadline, numpy was really not
        # collaborating, and eventually...
        MW_bulge_mean = []
        MW_halo_mean = []
        MW_disk_mean = []

        M31_bulge_mean = []
        M31_halo_mean = []
        M31_disk_mean = []

        M31_bulge_stdev = []
        M31_halo_stdev = []
        M31_disk_stdev = []

        MW_bulge_stdev = []
        MW_halo_stdev = []
        MW_disk_stdev = []

        for i in range(0, 800, 5):
            print(f"Snapsnot: {i}")
            MW_arr = np.zeros((3, 2))
            M31_arr = np.zeros((3, 2))
            try:
                MW = GalaxyPos(i, "MW")
                M31 = GalaxyPos(i, "M31")

                for idx, ptype in enumerate(ptypelist):
                    x, y, z = MW.position(ptype, COM=True)
                    x2, y2, z2 = M31.position(ptype, COM=True)

                    r = np.sqrt(x**2 + y**2 + z**2)
                    r2 = np.sqrt(x2**2 + y2**2 + z2**2)

                    MW_arr[idx][0] += np.mean(r)
                    M31_arr[idx][0] += np.mean(r2)
                    MW_arr[idx][1] += np.std(r)
                    M31_arr[idx][1] += np.std(r2)
            except:
                print("bad snapshot")

            MW_bulge_mean.append(MW_arr[2][0])
            M31_bulge_mean.append(M31_arr[2][0])

            MW_disk_mean.append(MW_arr[1][0])
            M31_disk_mean.append(M31_arr[1][0])
                
            MW_halo_mean.append(MW_arr[0][0])
            M31_halo_mean.append(M31_arr[0][0])
            
            MW_bulge_stdev.append(MW_arr[2][1])
            M31_bulge_stdev.append(M31_arr[2][1])

            MW_disk_stdev.append(MW_arr[1][1])
            M31_disk_stdev.append(M31_arr[1][1])
                
            MW_halo_stdev.append(MW_arr[0][1])
            M31_halo_stdev.append(M31_arr[0][1])
        
            ## truly cursed
            np.savetxt("meanstdev.txt", (MW_bulge_mean, MW_halo_mean, MW_disk_mean, M31_bulge_mean,    M31_halo_mean, M31_disk_mean, M31_bulge_stdev, M31_halo_stdev, M31_disk_stdev,     MW_bulge_stdev, MW_halo_stdev,     MW_disk_stdev))

        MW_bulge_mean, MW_halo_mean, MW_disk_mean, M31_bulge_mean,    M31_halo_mean, M31_disk_mean, M31_bulge_stdev, M31_halo_stdev, M31_disk_stdev,     MW_bulge_stdev, MW_halo_stdev,     MW_disk_stdev = np.loadtxt("meanstdev.txt")

        fig, ax = plt.subplots(2, ncols=len(ptypelist), figsize=[30, 15], sharex=True) #squeeze=False) #subplot_kw=dict(projection='3d'))
        j = np.array([i for i in range(0, 800, 5)])

        plt.rcParams.update({"font.size": 22})

        #print(j)
        #print(MW_halo_mean)

        ## milky way stuff
        ax[0, 0].plot(j, MW_halo_mean, label="Halo Mean")
        ax[0, 0].plot(j, MW_halo_stdev, label="Halo Stdev")
        ax[0, 0].legend(loc='center right')

        ax[0, 1].plot(j, MW_disk_mean, label="Disk Mean")
        ax[0, 1].plot(j, MW_disk_stdev, label="Disk Stdev")
        ax[0, 1].legend(loc='center right')

        ax[0, 2].plot(j, MW_bulge_mean, label="Bulge Mean")
        ax[0, 2].plot(j, MW_bulge_stdev, label="Bulge Stdev")
        ax[0, 2].legend(loc='center right')

        ## m31 stuff
        ax[1, 0].plot(j, M31_halo_mean,  label="Halo Mean")
        ax[1, 0].plot(j, M31_halo_stdev, label="Halo Stdev")
        ax[1, 0].legend(loc='center right')
                                       
        ax[1, 1].plot(j, M31_disk_mean,  label="Disk Mean")
        ax[1, 1].plot(j, M31_disk_stdev,  label="Disk Stdev")
        ax[1, 1].legend(loc='center right')
                                       
        ax[1, 2].plot(j, M31_bulge_mean, label="Bulge Mean")
        ax[1, 2].plot(j, M31_bulge_stdev, label="Bulge Stdev")
        ax[1, 2].legend(loc='center right')

        fig.supxlabel("Snapshot Number")
        fig.supylabel("Distance [kpc]")
        ax[0, 0].set_ylabel("Milky Way")
        ax[1, 0].set_ylabel("M31")
        fig.suptitle("MW and M31 radii Mean/Stdev")

        fig.savefig("meanstdev.png")


    ## showing the change in radii of MW and M31
    do_flag = False
    if do_flag:
        for i in range(0, 800, 5):
            try:
                fig = plt.figure(figsize=[30, 15])
                gs = GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1])

                MW = GalaxyPos(i, "MW")
                ax_gal = fig.add_subplot(gs[:, 1])
                mwx, mwy, mwz = MW.position(2, COM=True) ## only plotting disk
                hist = ax_gal.hist2d(mwx,mwy, bins=200, norm=LogNorm(), cmap='magma')
                ax_gal.set_xlabel("x [kpc]")
                ax_gal.set_ylabel("y [kpc]")
                #fig.colorbar(hist, label='Number  of  stars  per  bin')
                plt.colorbar(hist[3], label='Number  of  stars  per  bin')


                for idx, ptype in enumerate(ptypelist):
                    #idx = idx + 1 if idx > 0 else idx ## i hate matplotlib
                    ax = fig.add_subplot(gs[idx, 0])


                    ## change in radius over time
                    ## this is best illustrated through a histogram
                    ## showing the distribution of the radii of various 
                    ## ptypes. and how those distributions change over time.
                    ## we graph these against the calculated radius of the
                    ## galaxy, to get a sense for how many particles have left the galaxy
                    ## and also give us a sense of how far out they are
                    ## to find the radius, we look for the point where the bulge/disk levels off i.e. the inflection point of the total disk/bulge mass

                    mwx, mwy, mwz = MW.position(ptype, COM=True)
                    r = np.sqrt(mwx**2+mwy**2+mwz**2)    
                    ax.hist(r, bins=200)
                    ax.set_title(f"r distribution for MW, {ptypecheatsheet[ptype]}")
                    ax.set_xlabel("r [kpc]")
                    ax.set_ylabel("No. of particles at this distance")

                fig.savefig(f"output/histpos/MWRadii_{['000'+str(i)][0][-3:]}")
                plt.close()

                fig = plt.figure(figsize=[30, 15])
                gs = GridSpec(3, 2, width_ratios=[1, 1], height_ratios=[1, 1, 1])

                M31 = GalaxyPos(i, "M31")
                ax_gal = fig.add_subplot(gs[:, 1])
                m31x, m31y, m31z = MW.position(2, COM=True) ## only plotting disk
                hist = ax_gal.hist2d(m31x,m31y, bins=200, norm=LogNorm(), cmap='magma')
                ax_gal.set_xlabel("x [kpc]")
                ax_gal.set_ylabel("y [kpc]")
                plt.colorbar(hist[3], label='Number  of  stars  per  bin')

                for idx, ptype in enumerate(ptypelist):
                    ax = fig.add_subplot(gs[idx, 0])


                    ## change in radius over time
                    ## this is best illustrated through a histogram
                    ## showing the distribution of the radii of various 
                    ## ptypes. and how those distributions change over time.
                    ## we graph these against the calculated radius of the
                    ## galaxy, to get a sense for how many particles have left the galaxy
                    ## and also give us a sense of how far out they are
                    ## to find the radius, we look for the point where the bulge/disk levels off i.e. the inflection point of the total disk/bulge mass

                    m31x, m31y, m31z = M31.position(ptype, COM=True)
                    r = np.sqrt(m31x**2+m31y**2+m31z**2)    
                    ax.hist(r, bins=200)
                    ax.set_title(f"r distribution for M31 {ptypecheatsheet[ptype]}")
                    ax.set_xlabel("r [kpc]")
                    ax.set_ylabel("No. of particles at this distance")


                fig.savefig(f"output/histpos/M31Radii_{['000'+str(i)][0][-3:]}")
                plt.close()
            except:
                print(f"Problematic snapshot : {i}")


    ## this function would've dynamically detect the bounds of 
    ## the galaxies. I say "would've" because it doesn't work.
    do_flag = False
    if do_flag:

        for i in range(0, 800, 5):
                
            fig, ax = plt.subplots(figsize=[10, 5])
            for index, galaxy_name in enumerate(['M31', 'MW',], 0): 

                print(f"    {galaxy_name}")
                for idx, ptype in enumerate(ptypelist):
                    ## calculating r for each ptype
                    print("     Calculating r")
                    M = GalaxyPos(i, galaxy_name)
                    x, y, z = M.position(ptype, COM=True)
                    r = np.sqrt(x**2+y**2+z**2)    
                    ax.hist(r, bins=200)
                    ax.set_title(f"r distribution for {galaxy_name}, {ptypecheatsheet[ptype]}")

                    # calculating the inflection points (ito r in kpc)
                    ## (unfortunately, this routine does not work great
                        ## and was discarded)
                    print("     Calculating inflection points")
                    M = MassProfile(galaxy_name, i)

                    for j, pptype in {2:("Disk", "--"), 3:("Bulge", ":"),}.items():

                        print(f"        IFP for {pptype[0]}")
                        ## note: gaussian filter to smooth the function
                        particle_mass_distribution_smoothed = M.MassEnclosed(j, r) ## faussianfilter1d, 100

                        ## compute second derivative to find inflection point
                        pmds_second_derivative = np.gradient(np.gradient(particle_mass_distribution_smoothed))

                        ## compute inflection point
                        infls = np.where(np.diff(np.sign(pmds_second_derivative)))[0]
                        if len(infls) >= 1:
                            ax[index, idx].axvline(x=infls[0], label=f"Approximate Radius of {pptype[0]}")
                        ax[index, idx].legend(loc='upper right')

                fig.savefig(f"output/r_hist_{i}")

                
    ## function that generates all energies
    ## and plots the results
    do_flag = True
    if do_flag:
        
        for i in range(300, 500, 25):
            print(i)
                
            fig, ax = plt.subplots(4, 3, figsize=[30, 15])

            for idx, ptype in enumerate(ptypelist):

                EMW = GenerateEnergies(i, "MW", ptype)
                mw_to_mw, mw_to_m31 = EMW.check_if_bound()

                EM31 = GenerateEnergies(i, "M31", ptype)
                m31_to_mw, m31_to_m31 = EM31.check_if_bound()

                # to generate all potential files
                del EMW
                del EM31
                gc.collect()


                # Potentials w/ time.
                # this plot was not very useful...
                do_flag = False
                if do_flag: 
                    print("np.where")
                    print(len(np.where(mw_to_mw - mw_to_m31 > 0)[0]))
                    ax[idx].bar("Particles MW Bound", len(np.where(mw_to_mw - mw_to_m31 > 0)[0]), )
                    ax[idx].bar("MW Particles M31 Bound", len(np.where(mw_to_mw - mw_to_m31 < 0)[0]),)
                    ax[idx].bar("M31 Particles M31 Bound", len(np.where(m31_to_m31 - m31_to_mw > 0)[0]),)
                    ax[idx].bar("M31 Particles MW Bound", len(np.where(m31_to_m31 - m31_to_mw < 0)[0]),)

                    ax[idx].legend(loc='upper left')

                    print(f"mw {len(np.where(mw_to_mw<mw_to_m31)[0])}")
                    print(f"m31 {len(np.where(m31_to_m31<m31_to_mw)[0])}") ## interesting note: no unbound particles ever! which makes sense. they'll still be in the dm halo, and the halo kinematics don't change a whole lot during the collision.

                ## potentials energy histogram
                do_flag = False
                if do_flag:
                    nbins = 500
                    ax[0, idx].hist(mw_to_mw, bins=nbins)
                    ax[0, idx].set_title(f"MW {ptypecheatsheet[ptype]} Particles Energy wrt MW")
                    ax[3, 0].set_xlabel("Energy (Absolute Value) (m^3 / s^2 1e10Msun)")
                    ax[3, 0].set_ylabel("No. of particles at this Energy")
                    ax[1, idx].hist(mw_to_m31, bins=nbins)
                    ax[1, idx].set_title(f"MW {ptypecheatsheet[ptype]} Particles Energy wrt M31")
                    ax[2, idx].hist(m31_to_mw, bins=nbins)
                    ax[2, idx].set_title(f"M31 {ptypecheatsheet[ptype]} Particles Energy wrt MW")
                    ax[3, idx].hist(m31_to_m31, bins=nbins)
                    ax[3, idx].set_title(f"M31 {ptypecheatsheet[ptype]} Particles Energy wrt M31")
                

            fig, ax = plt.subplots(ncols=2, figsize=[30, 15])

            for idx, ptype in enumerate([2]):
                ## we only care about disk particles here
                
                EMW = GenerateEnergies(i, "MW", ptype)
                mw_to_mw, mw_to_m31 = EMW.check_if_bound()

                EM31 = GenerateEnergies(i, "M31", ptype)
                m31_to_mw, m31_to_m31 = EM31.check_if_bound()

                MW = GalaxyPos(i, "MW")
                M31 = GalaxyPos(i, "M31")

                ## position of all stuff
                
                ## transparent backgroud
                scat1 = ax[0].scatter(*MW.position(ptype=ptype),label="MW", color=[1., 0.5, 0.5])
                scat2 = ax[1].scatter(*M31.position(ptype=ptype),label="M31", color=[1., 0.5, 0.5])

                ## overlay
                scat1 = ax[0].scatter(*MW.position(ptype=ptype),label="MW", c=np.where(mw_to_mw-mw_to_m31<0, 'k', np.where(mw_to_mw-mw_to_m31>0, 'r', 'b'))) 
                scat2 = ax[1].scatter(*M31.position(ptype=ptype),label="M31", c=np.where(m31_to_m31-m31_to_mw<0, 'k', np.where(m31_to_m31-m31_to_mw>0, 'r', 'b')),)
                #ax[0, idx].legend()
                
                #scat1 = ax[1, idx].scatter(*MW.position(ptype=ptype),label="MW", c=(np.where(mw_to_mw>0)), cmap='magma') ## mw_to_mw - mw_to_m31???
                #ax[2, idx].set_title(f"Milky Way {ptypecheatsheet[ptype]}")
                #scat2 = ax[0, idx].scatter(*M31.position(ptype=ptype),label="M31", c=(np.where(m31_to_m31>0)), cmap='magma')
                #ax[2, idx].set_title(f"M31 {ptypecheatsheet[ptype]}")


            fig.savefig(f"output/potentialbar/{i}")



    do_flag = False
    if do_flag:
        for i in range(200, 300, 10):
            print(i)
                
            fig, ax = plt.subplots(ncols=3, figsize=[30, 15])

            for idx, ptype in enumerate(ptypelist):


                # plot positions of particles. to demonstrate mass transfer:
                # first plot MW on top of M31
                # then plot M31 on top of MW
                # (close encounter b/w 300 and 400)
                MW = GalaxyPos(i, "MW")
                M31 = GalaxyPos(i, "M31")
                mwpos = MW.position(ptype=ptype)
                m31pos = M31.position(ptype=ptype)

                ax1 = ax[idx].scatter(mwpos[0], mwpos[1], s=1, label="MW", c='b')#mw_to_mw-mw_to_m31, cmap='magma') 
                ax2 = ax[idx].scatter(m31pos[0], m31pos[1], s=1, label="M31", c='r')#m31_to_m31-m31_to_mw, cmap='magma')
                ax[idx].set_title(f"{ptypecheatsheet[ptype]}")
                ax[idx].set_xlabel("r [kpc]")
                ax[idx].set_ylabel("r [kpc]")

                ax[idx].legend(loc='upper left') ## THIS DOES NOT WORK!!! WHY???
                fig.suptitle(f"MW vs M31 (Snapshot {i})")


                ## 3d plotting : ax.view_init 

                #fig.colorbar(scat1, label='Boundness', ax=ax[0, idx])
                #fig.colorbar(scat2, label='Boundness', ax=ax[0, idx])

                #plt.title("Total Energy of Particles in Disks")

            fig.savefig(f"output/transfer/{i}")
        #plt.show()
        #plt.figlegend(["MW", "M31"])
        #fig.savefig(f"output/bound{i}")
    #ax.legend()
    #plt.show()

if __name__ == "__main__":
    main()

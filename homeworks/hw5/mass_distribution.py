from CenterOfMass import CenterOfMass
from readfile import Read
import astropy.units as u
import numpy as np
from astropy.constants import G
import matplotlib.pyplot as plt

G = G.to(u.kpc * u.km**2 / u.s**2 / u.M_sun)

class MassProfile:
    def __init__(self, galaxy, snapshot):
        """ Initializing the Mass Profile class.

        Inputs:
            galaxy: string
               Name of the galaxy we are examining.
            snapshot: int
                Snapshot number

        Outputs:
            None
        
        """
        self.filename = f'{galaxy}_{["000"+str(snapshot)][0][-3:]}.txt'
        self.gname = galaxy
        self.time, self.total, self.data = Read(self.filename)

        self.m = self.data['m']
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc

    def MassEnclosed(self, ptype, r):
        """ Function to return the mass of a given particle type within a given radius.

        Inputs:
            ptype: int
                Particle type. 1: Disk. 2: Disk. 3. Bulge.
            r: float or array of floats
                Distance from the center of the galaxy in kpc

        Outputs:
            mass: numpy array
                Array of masses in units of Msun
        """
        if type(r) == float:
            r = [r]
        COM = CenterOfMass(self.filename, ptype)
        COM_P = COM.COM_P(0.1)
        
        index = np.where(self.data['type'] == ptype)
        
        masses = []
        r2 = np.sqrt((self.x[index] - COM_P[0])**2+(self.y[index] - COM_P[1])**2+(self.z[index] - COM_P[2])**2)
        
        for i in r:
            index2 = np.where(r2.value <= i)
            masses.append(np.sum(self.m[index][index2]))

        return np.array(masses) * 1e10 * u.M_sun
    
    def MassEnclosedTotal(self, r,):
        """ Returns the total mass enclosed within a given radius.

        Inputs:
            r: float or array of floats
                Distance from the center of the galaxy in kpc.

        Outputs:
            mass: numpy array
                Array of masses in units of Msun

        """
        x = []
        for i in range(1, 4):
            if self.gname == "M33":
                if i == 3:
                    break
            if x == []:
                x = self.MassEnclosed(i, r)
            else:
                x += self.MassEnclosed(i, r)
        return np.array(x) * u.M_sun

    def HernquistMass(self, r, a, Mhalo):
        """ Function that defines Hernquist (1990) mass profile.
    
    Inputs:
        r: astropy quantity
            galactocentric distance in kpc
        
        a: astropy quantity
            scale radius of the Hernquist profile in kpc
        
        m_halo: float
            total halo mass in units of 1e12 M_sun
    
    Output:
        mass: astropy quantity
            total mass within the given radius in M_sun
        """
        return (Mhalo*r**2)/((a+r)**2)

    def CircularVelocity(self, ptype, r):
        """ Function to calculate the circular velocity of a given particle type at a given radius.
        Inputs:
            ptype: int
                Particle type. 1: Disk. 2: Disk. 3. Bulge.
            r: float or array of floats
                Distance from the center of the galaxy in kpc

        Outputs:
            velocities: numpy array
                Array of velocities in units of km/s

        """
        return np.round(np.sqrt(G*self.MassEnclosed(ptype, r)/r), 2)

    def TotalCircularVelocity(self, r):
        """  Function to calculate the circular velocity for all particles at a given radius.
        Inputs:
            r: float or array of floats
                Distance from the center of the galaxy in kpc

        Outputs:
            velocities: numpy array
                Array of velocities in units of km/s

        """
        return np.sqrt(G*self.MassEnclosedTotal(r)/r)

    def HernquistVCirc(self, r, a, Mhalo):
        """Function that returns velocities inside a given radius following a Hernquist (1990) mass profile.
    
    Inputs:
        r: astropy quantity
            galactocentric distance in kpc
        
        a: astropy quantity
            scale radius of the Hernquist profile in kpc
        
        m_halo: float
            total halo mass in units of 1e12 M_sun
    
    Output:
        velocity: astropy quantity
            velocity of particles in units of km/s
        """
        return np.round(np.sqrt(G*self.HernquistMass(r, a, Mhalo)/r), 2)

def main():
    fig, ax = plt.subplots(2, 3, sharey="row", sharex=True, figsize=(14, 6))

    ## iterate over the different galaxies
    for index, galaxy_name in enumerate(['M33', 'MW','M31',], 0): 
        
        ## Mass Calculation 
        print(f"Calculating for {galaxy_name}...")

        ## finding the rotation curves
        M = MassProfile(galaxy_name, 0)
        r = np.arange(0.1, 30, 0.1)
    
        for i, ptype in {1: ("Halo", "dashdot"), 2:("Disk", "--"), 3:("Bulge", ":"),}.items():
            if M.gname == "M33":
                if i == 3:
                    break

            print(f"	Mass for {M.gname}. Ptype: {ptype[0]}")
            ax[0, index].semilogy(r, M.MassEnclosed(i, r), linestyle=ptype[1], linewidth=3, 
                   label=ptype[0])
            print(f"	Velocity for {M.gname}. Ptype: {ptype[0]}")
            ax[1, index].semilogy(r, M.CircularVelocity(i, r), linestyle=ptype[1], linewidth=3, 
                   label=ptype[0])
        
        ax[1, index].semilogy(r, M.TotalCircularVelocity(r), linestyle=(5, (10, 3)), linewidth=3, 
                   label="Total")
        ax[0, index].semilogy(r, M.MassEnclosedTotal(r), linestyle=(5, (10, 3)), linewidth=3, 
                   label="Total")
        
        # setting up Hernquist 
        actual_mass = M.MassEnclosed(1, r)
        rmse_mass = None
        a_final_mass = None
        final_mass = []
        
        actual_velocity = M.CircularVelocity(1, r)
        rmse_velocity = None
        a_final_velocity = None
        final_velocity = []

        ## fitting the Hernquist Profile
        print(f"	Hernquist mass + velocity for {M.gname}.")
        for a in range(10):
            ## find masses and velocities predicted by the Hernquist profile
            full_halo_mass = M.MassEnclosed(1, 100.0)
            predicted_mass = M.HernquistMass(r, a, full_halo_mass)
            predicted_velocity = M.HernquistVCirc(r, a, full_halo_mass)

            ## normalized RMSE
            rmse_new_mass = np.linalg.norm(predicted_mass - actual_mass) / np.sqrt(len(actual_mass)) / np.mean(actual_mass) 
            rmse_new_velocity = np.linalg.norm(predicted_velocity - actual_velocity) / np.sqrt(len(actual_velocity)) / np.mean(actual_velocity) 
            
            ## find min rmse separately for velocity and mass
            if rmse_mass == None or rmse_new_mass < rmse_mass:
                rmse_mass = rmse_new_mass
                a_final_mass = a
                final_mass = predicted_mass
        
            if rmse_velocity == None or rmse_new_velocity < rmse_velocity:
                rmse_velocity = rmse_new_velocity
                a_final_velocity = a
                final_velocity = predicted_velocity
        

        ax[0, index].semilogy(r, final_mass,  color="fuchsia", linewidth=3, 
                   label=f"Hernquist.")

        ax[0, index].text(0.05, 0.95, f"{galaxy_name}", transform=ax[0, index].transAxes, size=11)
        ax[0, index].text(0.05, 0.90, f"Normalized RMSE: {rmse_mass.value:2e}.", transform=ax[0, index].transAxes, size=11)
        ax[0, index].text(0.05, 0.85, f"$a_{{final}}={a_final_mass}$", transform=ax[0, index].transAxes, size=11)
        
        
        ax[1, index].semilogy(r, final_velocity,  color="fuchsia", linewidth=3, 
                   label=f"Hernquist.")
        ax[1, index].text(0.15, 0.05, f"{galaxy_name}", transform=ax[1, index].transAxes, size=11)
        ax[1, index].text(0.15, 0.1, f"Normalized RMSE: {rmse_velocity.value:2e}.", transform=ax[1, index].transAxes, size=11)
        ax[1, index].text(0.15, 0.15, f"$a_{{final}}={a_final_velocity}$", transform=ax[1, index].transAxes, size=11)
    

    fig.suptitle("Mass Distributions and Velocity Curves in M31, M33, and the MW") 
    
    ax[1, 0].set_xlabel("R [kpc]")
    ax[0, 0].set_ylabel(r"Mass [$M_{sun}$]")
    ax[1, 0].set_ylabel("V [km/s]")
    
    plt.plot()
    ax[0, 2].legend(loc="lower right")
    ax[1, 2].legend(loc="lower right")
    
    ## special legends for M33
    ax[0, 0].legend(loc="lower right")
    ax[1, 0].legend(loc="upper right")

    plt.show()

if __name__ == "__main__":
    main()

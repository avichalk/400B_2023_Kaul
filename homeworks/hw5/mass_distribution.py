from CenterOfMass import CenterOfMass
from readfile import Read
import astropy.units as u
import numpy as np
from astropy.constants import G
import matplotlib.pyplot as plt


G = G.to(u.kpc * u.km**2 / u.s**2 / u.M_sun)


class MassProfile:
    def __init__(self, galaxy, snapshot):
        """

        """
        self.filename = f'{galaxy}_{["000"+str(snapshot)][0][-3:]}.txt'
        self.gname = galaxy
        self.time, self.total, self.data = Read(self.filename)

        self.m = self.data['m']
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc

    def MassEnclosed(self, ptype, r):
        """

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
        """

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
        """

        """
        return np.round(np.sqrt(G*self.MassEnclosed(ptype, r)/r), 2)

    def TotalCircularVelocity(self, r):
        """

        """
        return np.sqrt(G*self.MassEnclosedTotal(r)/r)

    def HernquistVCirc(self, r, a, Mhalo):
        """

        """
        return np.round(np.sqrt(G*self.HernquistMass(r, a, Mhalo)/r), 2)

def main():
    fig, ax = plt.subplots(2, 3, sharey="row", sharex=True, figsize=(14, 6))
    for index, galaxy_name in enumerate(['M33', 'MW','M31',], 0): 
        """
        Mass Calculation 
        """
        
        ## finding the rotation curves
        M = MassProfile(galaxy_name, 0)
        r = np.arange(0.1, 30, 0.1)
    
        for i, ptype in {2:("Disk", "--"), 3:("Bulge", ":"),}.items(): # 1:("Halo", "dashdot")}.items():
            if M.gname == "M33":
                if i == 3:
                    break
                #print(M.MassEnclosed(i, r)) #
            ax[0, index].semilogy(r, M.MassEnclosed(i, r), linestyle=ptype[1], linewidth=3, 
                   label=ptype[0])
        ax[0, index].semilogy(r, M.MassEnclosedTotal(r), linestyle=(5, (10, 3)), linewidth=3, 
                   label="Total Mass")
         
        ## finding the Hernquist Mass
        actual = M.MassEnclosed(1, r)
        rmse = None
        a_final = None
        final = []
        for a in range(100):
            predicted = M.HernquistMass(r, a, actual[-1])
            rmse_new = np.linalg.norm(predicted - actual) / np.sqrt(len(actual)) / np.mean(actual) ## normalized RMSE, though I wouldn't be surprised if dividing by the mean is incorrect. Either way, we find the minimum.
            if rmse == None or rmse_new < rmse:
                rmse = rmse_new
                a_final = a
                final = predicted
        
        ax[0, index].semilogy(r, actual, linestyle='dashdot', linewidth=3, 
                 label="Halo Mass")
        ax[0, index].semilogy(r, final,  color="fuchsia", linewidth=3, 
                   label=f"Hernquist Mass.")
        ax[0, index].text(0.1, 0.95, f"{galaxy_name}", transform=ax[0, index].transAxes, size=11)
        ax[0, index].text(0.1, 0.92, f"Normalized RMSE: {rmse.value:2e}.", transform=ax[0, index].transAxes, size=11)
        ax[0, index].text(0.1, 0.89, f"$a_{{final}}={a_final}$", transform=ax[0, index].transAxes, size=11)
        
        
        """
        Velocity Calculation
        """ 
        
        M = MassProfile(galaxy_name, 0)
        r = np.arange(0.1, 30, 0.1)
        
        ## finding the rotation curves
        for i, ptype in {3:("Bulge", ":"), 2:("Disk", "--"), 1:("Halo", "dashdot")}.items():
            if M.gname == "M33":
                if i == 3:
                    break
            ax[1, index].semilogy(r, M.CircularVelocity(i, r), linestyle=ptype[1], linewidth=3, 
                   label=ptype[0])
        ax[1, index].semilogy(r, M.TotalCircularVelocity(r), linestyle=(5, (10, 3)), linewidth=3, 
                   label="Total Mass")
      
        ## Fitting a Hernquist Profile 
        actual = M.CircularVelocity(1, r)
        rmse = None
        a_final = None
        final = []
        for a in range(100):
            predicted = M.HernquistVCirc(r, a, actual[-1])
            rmse_new = np.linalg.norm(predicted.value - actual.value) / np.sqrt(len(actual)) / np.mean(actual) ## normalized RMSE, though I wouldn't be surprised if dividing by the mean is incorrect. Either way, we find the minimum.
            if rmse == None or rmse_new < rmse:
                rmse = rmse_new
                a_final = a
                final = predicted
        
        # ax[1, index].semilogy(r, actual, linestyle='dashdot', linewidth=3, 
        #         label="Halo Mass")
        ax[1, index].semilogy(r, final,  color="fuchsia", linewidth=3, 
                   label=f"Hernquist Mass.")
        ax[1, index].text(0.2, 0.55, f"{galaxy_name}", transform=ax[1, index].transAxes, size=11)
        ax[1, index].text(0.2, 0.52, f"Normalized RMSE: {rmse.value:2e}.", transform=ax[1, index].transAxes, size=11)
        ax[1, index].text(0.2, 0.49, f"$a_{{final}}={a_final}$", transform=ax[1, index].transAxes, size=11)
    

    fig.suptitle("Mass Distributions and Velocity Curves in M31, M33, and the MW") 
    
    ax[1, 0].set_xlabel("R [kpc]")
    ax[0, 0].set_ylabel(r"Mass [$M_{sun}$]")
    ax[1, 0].set_ylabel("V [km/s]")
    
    plt.plot()
    ax[0, 2].legend(loc="lower right")
    ax[1, 2].legend(loc="center right")
    plt.show()

if __name__ == "__main__":
    main()

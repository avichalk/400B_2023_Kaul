
# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# Make edits where instructed - look for "****", which indicates where you need to 
# add code. 




# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from galaxy_mass import component_mass

# # M33AnalyticOrbit


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): # **** add inputs
        """ **** ADD COMMENTS """

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        ### **** store the output file name
        self.filename = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        M33_COM = CenterOfMass("M33_000.txt", 2)

        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_P = M33_COM.COM_P(0.1)
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_COM_V = M33_COM.COM_V(*M33_COM_P).value
        
        M31_COM = CenterOfMass("M31_000.txt", 2)
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)

        M31_COM_P = M31_COM.COM_P(0.1)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_COM_V = M31_COM.COM_V(*M31_COM_P).value
        
        M31_COM_P, M33_COM_P = M31_COM_P.value, M33_COM_P.value
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        self.r0 = M33_COM_P - M31_COM_P 
        self.v0 = M33_COM_V - M31_COM_V

        ### get the mass of each component in M31 
        ### disk (ptype 2)
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 ## kpc

        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = component_mass("M31_000.txt", 2) * 1e12

        ### bulge (ptype 3)
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1 ## kpc
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = component_mass("M31_000.txt", 3)

        # Halo (ptype 1)
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 ## kpc. from Assignment 4 scale length
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = component_mass("M31_000.txt", 1)
    
    
    def HernquistAccel(self, M, rhalo, r): # it is easiest if you take as an input the position VECTOR 
        """ **** ADD COMMENTS """
        
        ### **** Store the magnitude of the position vector
        rmag = np.sqrt(r.dot(r))
        
        ### *** Store the Acceleration
        Hern = - self.G * M / rmag * 1 / (rhalo + rmag) ** 2 * r
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, r_d, r):# it is easiest if you take as an input a position VECTOR  r 
        """ **** ADD COMMENTS """

        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        
        R = np.sqrt(r[0]**2 + r[1]**2)
        B = self.rdisk + np.sqrt(r[2]**2 + self.rdisk)
        a = - self.G * M * r / (R**2 + B**2) ** 1.5 \
                * np.array([1, 1, B/np.sqrt(r[2]**2 + self.rdisk**2)])
        
       
        return a
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, r): # input should include the position vector, r
        """ **** ADD COMMENTS """

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
            
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        a_bulge_halo = self.HernquistAccel((self.Mhalo + self.Mbulge), self.rhalo, r)
        a_disk = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, r)

        return np.sqrt(a_bulge_halo.dot(a_disk))
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """ **** ADD COMMENTS """
        
        # predict the position at the next half timestep
        rhalf = r + v * dt/2
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + self.M31Accel(rhalf) * dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = rhalf + vnew * dt/2
        
        return rnew, vnew # **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ **** ADD COMMENTS """
        
        # initialize the time to the input starting time
        t = t0
        
        # initialize an empty array of size :  rows int(tmax/dt)+2  , columns 7
        orbit = np.zeros((int(tmax/dt)+2, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while t<tmax:  # as long as t has not exceeded the maximal time 
            t += dt
            # **** advance the time by one timestep, dt
           
            # **** store the new time in the first column of the ith row
            orbit[i, 1] = t
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            
            pos, vel = self.LeapFrog(dt, orbit[i-1, 1:4], orbit[i-1, 4:])
    
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            orbit[i, 1:4] = pos
            orbit[i, 4:] = vel
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1 

            print(f"{t} of {tmax}")
        
        
        # write the data to a file
        np.savetxt(self.filename, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function

def main():
    #M33 = M33AnalyticOrbit("outfile.txt", )
    #M33.OrbitIntegration(0, 0.1, 10)


    # Read in the data files for the orbits of each galaxy that you just created
    # headers:  t, x, y, z, vx, vy, vz
    # using np.genfromtxt
    
    M31 = np.genfromtxt("output/LowRes/Orbit-M31.txt")
    M31_M33 = np.genfromtxt("output/LowRes/Orbit-M33.txt")

    # function to compute the magnitude of the difference between two vectors 
    # You can use this function to return both the relative position and relative velocity for two 
    # galaxies over the entire orbit  
    
    time = M31[:, 0]

    ## MW vs M31
    
    # Determine the magnitude of the relative position and velocities 
    
    # of MW and M31
    M31_Position = np.array([M31[:, 1], M31[:, 2], M31[:, 3]],)
    M31_Velocity = np.array([M31[:, 4], M31[:, 5], M31[:, 6]],)

    # of M33 and M31
    M33_M31_Separation = np.array([M31_M33[:, 1], M31_M33[:, 2], M31_M33[:, 3]],)
    M33_M31_Velocity = np.array([M31_M33[:, 4], M31_M33[:, 5], M31_M33[:, 6]],)

    
    fig, ax = plt.subplots(2, 2, sharey="row", sharex=True, figsize=(8, 4))
    
    # Plot the Orbit of the galaxies 
    #################################
    ax[0].plot(time, M31_Position, "MW vs M31 (Separation)")

    ax[0].plot(time, M33_M31_Separation)

    # Plot the orbital velocities of the galaxies
    #################################
    ax[1].plot(time, M31_Velocity)

    ax[1].plot(time, M33_M31_Velocity)

    ax[1].set_xlabel("Time [Gyr]")
    ax[0].set_ylabel("Separation [kpc]")
    ax[1].set_ylabel("V [km/s]")

    plt.savefig("fig.png")


if __name__ == "__main__":
    main()

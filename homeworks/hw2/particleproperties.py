import numpy as np
import astropy.units as u
from readfile import Read

def main():
    m, d, v = particleinfo("MW_000.txt", 'disk', 100)
    print(f'Mass: {m}, Distance: {d}, Velocity: {v}')

def particleinfo(name, ptype, no):
    """
    Takes the filename, data type, and particle number, and returns given properties of the relavant particle.
    Inputs: filename, type of particle (as a string), the particle no. you want (not in index format)
    """
    no -= 1 ## for array indexing purposes
    t, pno, data = Read(name) ## read data from file
    
    ## finding indices for diff. particle types.
    ## there is absolutely a better way to do this,
    ## but i can't figure it out
    darkin = np.where(data['type'] == 1)
    diskin = np.where(data['type'] == 2)
    bulgin = np.where(data['type'] == 3)
    
    ## particle properties depending on particle type
    ## particle type (int), mass (float, 10**10 M_sun), x, y, z (float, kpc), vx, vy, vz (float, km/s)
    if ptype=="dark":
        typ, m, x, y, z, vx, vy, vz = data[darkin][no]
    if ptype=="disk":
        typ, m, x, y, z, vx, vy, vz = data[diskin][no]
    if ptype=='bulge':
        typ, m, x, y, z, vx, vy, vz = data[bulgin][no]
    
    ## converting components into single 3d magnitude
    ## and then rounding values and assigning units
    ## the smart thing to do would be to assign them all
    ## in one block. but instead I did them all in different places
    d = np.around((components_to_3d(x, y, z)* u.kpc).to(u.lyr), 3) ## convert to one mag, assign unit kpc, convert to lyr, and round to 3 decimal places
    v = np.around(components_to_3d(vx, vy, vz), 3) * u.km / u.s ## convert to one mag, round to 3 decimal places, assign unit km/s
    return np.around(m * 10**10, 3) * u.M_sun, d, v ## assign unit 10**10 Msun to mass, round to 3 decimal places

def components_to_3d(x, y, z):
    ## Converts all 2d separated components into one 3d magnitude vector
    return np.sqrt(x**2 + y**2 + z**2)

if __name__ == "__main__":
    main()

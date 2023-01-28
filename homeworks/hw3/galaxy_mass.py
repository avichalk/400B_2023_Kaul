from readfile import Read
import numpy as np
import astropy.units as u
import pandas as pd
from tabulate import tabulate

def component_mass(filename, ptype,):
    t, pno, data = Read(filename) ## read data from file
    
    ## finding indices for diff. particle types.
    ## there is absolutely a better way to do this,
    ## but i can't figure it out
    darkin = np.where(data['type'] == 1)
    diskin = np.where(data['type'] == 2)
    bulgin = np.where(data['type'] == 3)
    
    ## particle properties depending on particle type
    ## particle type (int), mass (float, 10**10 M_sun), x, y, z (float, kpc), vx, vy, vz (float, km/s)
    mass = 0
    if ptype=="dark":
        for i in data[darkin]:
            mass += i[1] ## sums the masses of all particles belonging to this field
    if ptype=="disk":
        for i in data[diskin]:
            mass += i[1]
    if ptype=='bulge':
        for i in data[bulgin]:
            mass += i[1]

    return np.round(mass / 1e2, 3) ## rounds masses, as they are currently in 10^10Msun we return in units of 10^12 Msun by dividing by 1e2

def main():
    output = []
    local_group_mass = 0
    local_group_stellar_mass = 0
    for i in ["MW_000.txt", "M31_000.txt", "M33_000.txt"]:
        dark = component_mass(i, "dark")
        disk = component_mass(i, "disk")
        bulge = 0
        try: 
            bulge = component_mass(i, "bulge")
        except:
            pass
        total = dark + disk + bulge
        local_group_mass += total
        local_group_stellar_mass += disk + bulge

        output.append([i.split("_")[0], dark, disk, bulge if bulge!=0 else "N/A", np.round(total, 3), np.round((bulge+disk)/total, 3)])
   
    ## this separate table is one of the worst solutions to a fomatting problem I've ever come up with
    ## i feel dirty just looking at it.
    out2 = []
    out2.append([f"Total Mass of the Local Group [$10^{{12}} M_\odot$]: {np.round(local_group_mass, 3)}"])
    out2.append([f"$f_{{bar}}$ of the Local Group: {np.round(local_group_stellar_mass/local_group_mass, 3)}"])
    
    ## define headers (with tex formatting) and save to csv file
    headers = ["Galaxy Name", "Halo Mass [$10^{12} M_\odot$]", "Disk Mass [$10^{12} M_\odot$]", "Bulge Mass [$10^{12} M_\odot$]", "Total Mass [$10^{12} M_\odot$]", "$f_{bar}$"]
    pd.DataFrame(output, columns=headers).to_csv("output.txt")
    
    ## why spend two minutes copying over the contents of the csv file and formatting them into a tex document when you can spend half an hour automating the process?
    with open("tex/out.tex", "w") as f:
        f.writelines("\\documentclass{article} \n \\usepackage[a4paper, total={8in, 8in}]{geometry} \n \\begin{document} \n")
        f.writelines(tabulate(output, headers=headers, tablefmt='latex_raw'))
        f.writelines('\n\n')
        f.writelines(tabulate(out2, tablefmt='latex_raw'))
        f.writelines(f"\n\n 1. Total Mass of MW vs M31: Ratio of {np.round(output[0][4]/output[1][4], 3)}. {headers[np.argmax([output[0][i]/output[1][i] for i in range(1, 3)])+1].split('[')[0]}dominates.")
        f.writelines(f"\n\n 2. Stellar Mass of MW vs M31: Ratio of {np.round((output[0][2]+output[0][3])/(output[1][2]+output[1][3]), 3)}. As {['MW' if (output[0][2]+output[0][3])/(output[1][2]+output[1][3])>1 else 'M31'][0]} contains more stellar mass, which is the only type of mass that directly contributes to the luminosity of a galaxy, we expect it to be more luminous.")
        f.writelines(f"\n\n 3. Dark Matter Content of MW vs M31: Ratio of {np.round(output[0][1]/output[1][1])}. Yes - despite the stark difference in stellar matter, they have almost exactly the same dark matter content!")
        f.writelines(f"\n\n 4. Baryonic Fraction $f_{{bar}}$ for each galaxy is {', '.join(str(x) for x in [f'{output[i][0]}:{output[i][5]}' for i in range(0, 3)])}. Our baryonic fractions are universally lower than the provided value of 16\%. {', '.join(str(x) for x in [f'{output[i][0]}:{output[i][5]/.16}' for i in range(0, 3)])}. This discrepancy might come down to [TODO!]")
        f.writelines("\n \\end{document}")

if __name__ == "__main__":
    main()

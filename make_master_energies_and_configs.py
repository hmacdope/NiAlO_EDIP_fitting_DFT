import numpy as np
import pandas as pd
import os

def grab_final_config():
    with open("final_state.xyz", 'r') as f:
        lines = f.readlines()
        return lines

def add_tag(lines, dir):
    energy = float(lines[1].split()[1])
    lines[1] = lines[1].replace('\n','') + ' config: ' + dir + '\n'
    return lines, energy
    
def extract_eam_energy():
    with open("log.lammps", 'r') as f:
        lines = f.readlines()
        min_ene_line = 0
        for i, line in enumerate(lines):
            if '  Energy initial, next-to-last, final =' in line:
                min_ene_line = i + 1
                break
        if not min_ene_line:
            raise Exception("minimization_energy not found")
        eline = lines[min_ene_line]
        eval = float(eline.split()[0])
    return eval



dirs = [d for d in os.listdir('.') if os.path.isdir(d) and d != '.git' ]
dirs = sorted(dirs)

all_configs = []
energies  = {}
dft_energies = []
for dir in dirs:
    os.chdir(dir)
    config_lines = grab_final_config()
    tagged_lines, dft_energy = add_tag(config_lines, dir)
    dft_energies.append(dft_energy)
    all_configs += tagged_lines
    os.chdir('./EAM_energy')
    try:
        energies[dir] = extract_eam_energy()
    except:
        print(f"dir {dir} is likely an O only dir")
        energies[dir] = 0

    os.chdir('../../')

dirs = []
e_eam = []

for k,v in energies.items():
    dirs.append(k)
    e_eam.append(v)

df = pd.DataFrame({"Species": dirs, "EAM_Energies": e_eam, "DFT_Energies": dft_energies})
df["E_diff"] = df["DFT_Energies"] - df["EAM_Energies"]
df = df.set_index("Species")
print(df)

df.to_csv("Combined_data.csv")
df["EAM_Energies"].to_csv("EAM_energies.csv")
df["DFT_Energies"].to_csv("DFT_energies.csv")

with open("master_configs.xyz", 'w') as f:
    for line in all_configs:
        f.write(line)
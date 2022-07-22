import os 

flist = os.listdir()
file = [ f for f in flist if "no_ox" in f]
assert(len(file) == 1) 
file = file[0]

metals = ""
if "al" in file.lower():
    metals += "Al "
if "ni" in file.lower():
    metals += "Ni "



inp = f"""
################################################################################
# LAMMPS input file                                                            #
#                                                                              #
# Prototype run for single point potential                                     #
#                                                                              #
# richard terrett, 2018                                                        #
#                                                                              #
################################################################################

# System init

units      metal
dimension  3
boundary   p p p
atom_style atomic
timestep   0.01

# read data
box tilt large
read_data {file} 

# pair potential

pair_style eam/alloy
pair_coeff * * Mishin-Ni-Al-2009.eam.alloy {metals} 

# reporting

thermo 100

# Single point (0 iteration opt)

min_style cg
minimize 1e-5 1e-5 0 0

# Thermalisation

velocity all create 300 1 mom yes rot no

dump   1 all atom 10 generic.dump
fix 1 all nve
run  10000

write_restart checkpoint.init
"""

with open("single_point.run", 'w') as f:
    f.write(inp)

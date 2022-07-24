
# FITTING AN EDIP POTENTIAL FOR NiAl RNCs with adsorbed Oxygen

This directory contains DFT data and tools required to fit and EDIP potential
for NiAl reactive nanocomposites with adsorbed oxygen.

The corresponding potential used for The Ni-Al interaction is the 2009 Mishin Ni Al EAM/alloy potential

Each folder contains a DFT optimisation of a training system for fitting the potential.

Basic steps followed here are as follows:

1. Calculate DFT energy using VASP
2. Extract the **LAST** optimization step as a LAMMPS DATA file using PyIron/ASE in the notebook provided
3. Extract the **LAST** optimization step as above with the Oxygen stripped out
4. Use scripts to setup and run a LAMMPS calculation using the Mishin EAM potential to obtain EAM energy of the snapshot. Energy obtained by doing a 0 step minimzation
5. Use scripts to compile all the configurations and energies into the requisite files to then feed to the FORTRAN fitting routine 

Various scripts are provided in the folder to assist with this and should be run in this order

1. extract-all-states.sh grabs all the configurations from the optimization and writes XYZ file with box as extra atoms (dependent on a fortran routine extract_vasp.f90) 
2. extract_last.sh calls the python script extract_last_config.py to grab the last configuration from the above and write as an xyz file
3. extract_to_data.ipynb can read the VASP DFT calculations and spit out LAMMPS data files with and without oxygen. MUST BE RUN in a seperate folder above this one eg ../convert_format
4. create_EAM_folders.sh makes the EAM_energy subfolders within each DFT calculation directory and copies the DATA files created at step 2 there.
5. make_runfiles.sh calls the python script make_run.py in each EAM_energy subfolder to generate the LAMMPS input files
6. run-all-lammps.sh runs the EAM simulations in each subfolder
7. make_master_energies_and_configs.py extracts all the configurations as one big file and generates a bunch of CSVs with the energies.
8. fit!!! 



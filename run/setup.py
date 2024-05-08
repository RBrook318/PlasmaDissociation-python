#########################################################################################
#
#   Python Run script for setting up run paramters for dissociation after electron impact
#   Written by O Bramley                                          08.05.24
# 
#     
#
#########################################################################################

import os
import sys
import shutil
import subprocess
import reads_writes
import org_conv_create
import json

with open('inputs.json') as f:
    inputs=json.load(f)
# Runfolder name
Runfolder=inputs["setup"]["Runfolder"]


if os.path.exists("../"+Runfolder):
    value=input("File already exists do you want to delete it or continue with part two of setup? y/2/n\n")
    if(value=='y'):
        shutil.rmtree("../"+Runfolder)
        os.mkdir("../"+Runfolder)
        EXDIR="../"+Runfolder
        # Copies inputs.json file to execution folder
        shutil.copy2("inputs.json",EXDIR)
        # Copy run script to start folder
        shutil.copy2("run.py",EXDIR)
        # Copy restart script to start folder
        shutil.copy2("restart.py",EXDIR)
        # Name of geometry file
        Geometry=inputs["setup"]["Geometry_file"]
        atoms=["run"]["Atoms"]

        # Make Geom folder
        os.mkdir(EXDIR+"/Geom")
        # Makes folder for initial geometry
        os.mkdir(EXDIR+"/geom-"+Runfolder)
        geom_folder=EXDIR+"/geom-"+Runfolder

        # Read in geometry file
        coords=[]
        reads_writes.read_xyz(Geometry, atoms, coords)
        # Write the mass file file
        reads_writes.mass_write(EXDIR, atoms, coords)
        # Write Qchem input file
        reads_writes.write_opt_freq_input(geom_folder, atoms, coords)

        value=input("Run the geometry calculation as a job or locally? y/n\n")
        if(value=='y'):
            shutil.copy2("setup.py",EXDIR)
            # Calls subroutine that writes the script to run the qchem job
            reads_writes.run_script_write(geom_folder) 
            subprocess.call(["qsub", "geom.sh"],cwd=geom_folder) 
        else:
            subprocess.call(["qchem", "-save", "-nt", "8", "opt_freq.inp", "f.out", "wf"],cwd=geom_folder)
    elif(value=='2'):
        print("Continuing with part two of setup")
        val1=input("Is the runfolder called, "+ Runfolder +"? y/n\n")
        if(val1=='n'):
            val2=input("Type in the runfolder now \n")
            if os.path.exists("../"+val2):
                print("Folder exists proceeding")
                Runfolder=val2
            else: 
                sys.exit("Folder does not exist. Please check spelling or create folder")
        EXDIR="../"+Runfolder
        geom_folder=EXDIR+"/geom-"+Runfolder
        atoms=["run"]["Atoms"]
    else:
        sys.exit("Runfolder already exists. Change the Runfolder name or delte/move it")

mom_num=inputs["setup"]["repeats"]
T=inputs["run"]["Temp"]
#  Reads in f_out file
opt_geoms, modes, mode_cnt=reads_writes.read_f_out(geom_folder,atoms)
# Convert Geometries to bohr
opt_geoms=org_conv_create.convert_to_bohr(opt_geoms)
# Convert modes to form Create_geom.py uses
modes=org_conv_create.organise_modes(modes)
# Reads and multiplies by 1822.887
masses=reads_writes.read_mass()
# Creates the momenta
Px, Py, Pz = org_conv_create.create_geom(atoms,mode_cnt,T,modes,masses,mom_num)
# Writes momenta files to Geom folder
reads_writes.write_momentas(opt_geoms,(EXDIR+"/Geom"),Px,Py,Pz,atoms,mom_num)

# Qchem.py 15/11/2023

import os
import numpy as np 
import re
import math
from pyqchem import QchemInput, Structure
from pyqchem import get_output_from_qchem
np.set_printoptions(precision =30)

def file_contains_string(file_path, search_string):
    with open(file_path, "r") as file:
        for line in file:
            if search_string in line:
                return True
    return False

def create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS", Guess=True):
    
    # Filter indices based on dissociation flag
    active_indices = [i for i, flag in enumerate(molecule.dissociation_flags) if flag == 'NO']
    active_coords = [molecule.coordinates[i] for i in active_indices]
    active_symbols = [molecule.symbols[i] for i in active_indices]
    molecule = Structure(coordinates=active_coords, symbols=active_symbols, multiplicity=molecule.multiplicity)
   
    if spin_flip==0:
        qc_inp=QchemInput(molecule,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true'},
                        )       
    elif spin_flip==1:                
        qc_inp=QchemInput(molecule,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500},
                        # set_iter=500,
                        cis_n_roots=1,
                        cis_state_deriv=1
                        )
    if Guess:
       qc_inp.update_input({'scf_guess':'Read'})  
 
    return qc_inp
                      
def run_qchem(ncpu, molecule, n, nstates, spin_flip, Guess=True): 
    qc_inp=create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS", Guess=Guess)
    try:
        output = get_output_from_qchem(qc_inp,processors=ncpu)
    except:
        print('Using DIIS_GDM algorithm')
        # Retry with a different setup
        qc_inp.update_input({'scf_algorithm': 'DIIS_GDM', 'scf_guess': 'sad'})
        try:
            output = get_output_from_qchem(qc_inp,processors=ncpu)
        except:
            with open("ERROR", "w") as file:
                file.write("Error occurred during QChem job. Help.\n" + os.getcwd())
            exit()
    # Job completed successfully
    readqchem(output, molecule, n, nstates,spin_flip)
    # Append f.out content to f.all
    with open("f.all", "a") as f_all:
        f_all.write(output)
    return 

def readqchem(output, molecule, natoms, nst,spin_flip):
    reduced_natoms = sum(flag.lower() != 'yes' for flag in molecule.dissociation_flags)
    ndim = 3 * reduced_natoms
    if spin_flip==1:
        enum = output.find('Total energy for state  1:')
        scf_erg = float(output[enum: enum+100].split()[5])
        molecule.scf_energy[0]=scf_erg
        # molecule.scf_energy[1]=scf_erg
        l2t = ' Gradient of the state energy (including CIS Excitation Energy)'
    else:
        enum = output.find('Total energy in the final basis set')
        scf_erg = float(output[enum: enum+100].split()[8])
        molecule.scf_energy[0]=scf_erg
        # molecule.scf_energy[1]=scf_erg
        l2t = ' Gradient of SCF Energy'
    
    output_lines = output.split("\n")
    enum = output.find(l2t)
    output_lines = output[enum:-1].split("\n")
    lines_to_read = 4 * (math.ceil(reduced_natoms / 6)) +1
    forces = output_lines[1:lines_to_read]
    forces = [line.split() for line in forces]
    f = np.zeros(ndim,dtype = np.float64)
    strt = 0
    for i in range(int(len(forces)/4)):
        num=len(forces[i*4])
        for j in range(3):
            f[strt:strt+num] = forces[i*4+j+1][1:]
            strt = strt + num
    f = -f
    f = np.where(f == -0.0, 0.0, f)
    # Update the forces in the Molecule object
    molecule.update_forces(f)
    
   
# Qchem.py 15/11/2023

import os
import numpy as np 
import math
from pyqchem import QchemInput, Structure
from pyqchem import get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_frequencies import basic_frequencies

np.set_printoptions(precision =30)

# Electronic structure via QChem.
        
def file_contains_string(file_path, search_string):
    with open(file_path, "r") as file:
        for line in file:
            if search_string in line:
                return True
    return False

def create_qchem_input(molecule, spin_flip,scf_algorithm="DIIS", Guess=True):
    
    # Filter indices based on dissociation flag
    active_indices = [i for i, flag in enumerate(molecule.dissociation_flags) if flag == 'NO']
    active_coords = [molecule.coordinates[i] for i in active_indices]
    active_symbols = [molecule.symbols[i] for i in active_indices]
    coefficients = molecule.elecinfo
    molecule = Structure(coordinates=active_coords, symbols=active_symbols, multiplicity=molecule.multiplicity)
    if spin_flip==0:
        if Guess:             
            qc_inp=QchemInput(molecule,
                        scf_guess = coefficients,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true'}
                        # set_iter=500,
                        )
        else:
            qc_inp=QchemInput(molecule,
                        jobtype='force',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm=scf_algorithm,
                        extra_rem_keywords={'input_bohr':'true'}
                        # set_iter=500,
                        )  
    elif spin_flip==1:   
        if Guess:             
            qc_inp=QchemInput(molecule,
                        scf_guess = coefficients,
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
        else:
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
       
    return qc_inp
                      
def run_qchem(ncpu, molecule, n, nstates, spin_flip, Guess=True): 
    qc_inp=create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS", Guess=Guess)
    try:
        output, ee = get_output_from_qchem(qc_inp,processors=ncpu,return_electronic_structure=True)
        molecule.elecinfo=(ee['coefficients'])
    except:
        print('Using DIIS_GDM algorithm')
        # Retry with a different setup
        qc_inp=create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS_GDM", Guess=False)
        try:
            output, ee = get_output_from_qchem(qc_inp,processors=ncpu,return_electronic_structure=True)
            molecule.elecinfo=(ee['coefficients'])
        except:
            with open("ERROR", "w") as file:
                file.write("Error occurred during QChem job. Help.\n" + os.getcwd())
                file.write(output)
            exit()
    # Job completed successfully
    readqchem(output, molecule, n, nstates,spin_flip)
    # Append f.out content to f.all
    with open("f.all", "a") as f_all:
        f_all.write(output)
    return molecule


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
    for i in range(reduced_natoms):
        x = float(forces[1 + 4 * (i // 6)][i % 6 + 1])
        y = float(forces[2 + 4 * (i // 6)][i % 6 + 1])
        z = float(forces[3 + 4 * (i // 6)][i % 6 + 1])
        f[strt:strt+3] = [x, y, z]
        strt += 3
    f = -f
    f = np.where(f == -0.0, 0.0, f)
    # Update the forces in the Molecule object
    molecule.update_forces(f)

def initial_conditions(symbols,coords,cores):
    molecule = Structure(coordinates=coords, symbols=symbols, multiplicity=1)
    qc_inp = QchemInput(molecule,
                        jobtype='opt',
                        exchange='BHHLYP !50% HF +  50% Becke88 exchange',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm='diis',
                        extra_rem_keywords={'input_bohr':'true'}
                        )
    
    output = get_output_from_qchem(qc_inp,processors=cores)
    with open('optimisation.out', 'w') as f:
        f.write(output)
    parser_output = basic_optimization(output)
    opt_geoms=(parser_output['optimized_molecule'])
    # Find Normal modes
    qc_inp = QchemInput(opt_geoms,
                        jobtype='FREQ',
                        exchange='BHHLYP',
                        basis='6-31+G*',
                        extra_rem_keywords={'input_bohr':'true'}
                        )
    output = get_output_from_qchem(qc_inp,cores)
    with open('modes.out', 'w') as f:
        f.write(output)
    parser_output = basic_frequencies(output)
    return opt_geoms,parser_output
# -------------------------------------------------------------------------



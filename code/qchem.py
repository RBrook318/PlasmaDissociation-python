# Qchem.py 15/11/2023
import subprocess
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
                        scf_guess=Guess
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
                        scf_guess=Guess,
                        extra_rem_keywords={'input_bohr':'true','spin_flip':'true'},
                        set_iter=500
                        )
        return qc_inp
                      
    
    
    



def run_qchem(ncpu,file, molecule, n, nstates, spin_flip, Guess=True): 

    def submit_qchem_job(ncpu,file,output):
        subprocess.run(["qchem", "-save", "-nt", str(ncpu), file, output, "wf"])

    # Prepare f.inp file
    qc_inp=create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS", Guess=Guess)
    output = get_output_from_qchem(qc_inp,processors=ncpu)
   
   
    if "Total job time" in output:
        # Job completed successfully
        readqchem('f.out', molecule, n, nstates,spin_flip)
        # Append f.out content to f.all
        with open("f.out", "r") as f_out, open("f.all", "a") as f_all:
            f_all.write(f_out.read())
        pass   
    else:
        # Retry with a different setup
        qc_inp.update_input({'scf_algorithm': 'DIIS_GDM', 'scf_guess': 'false'})
        output = get_output_from_qchem(qc_inp,processors=ncpu)
        

        # Check for the second failure
        if "Total job time" in output:
            readqchem('f2.out', molecule, n, nstates,spin_flip)
            with open("f.out", "r") as f_out, open("f.all", "a") as f_all:
                f_all.write(f_out.read())
        else:
            with open("ERROR", "w") as file:
                file.write("Error occurred during QChem job. Help.\n" + os.getcwd())
           

    

def readqchem(output_file, molecule, natoms, nst,spin_flip):

    reduced_natoms = sum(flag.lower() != 'yes' for flag in molecule.dissociation_flags)
    ndim = 3 * reduced_natoms
    if spin_flip==1:
        l1t = ' Excited state   1: excitation energy (eV) ='
        l2t = ' Gradient of the state energy (including CIS Excitation Energy)'
    else:
        l1t = ' SCF   energy in the final basis set ='
        l2t = ' Gradient of SCF Energy'
    
    with open(output_file, 'r') as file:
        f = np.zeros(ndim,dtype = np.float64)

        molecule.update_forces(f)
        e = np.zeros(nst,dtype = np.float64)
        C = np.zeros(ndim)
        l1_found = False
        l2_found = False
        found_target = False
        for line in file:
            if found_target:
                # Read the line below the target line
                data_line = line.strip()
                match = re.search(r'-?\d+\.\d+', data_line)
                if match:
                    e[0] = float(match.group())
                    # Update the SCF energy in the Molecule object
                    molecule.update_scf_energy(e)

                else:
                    print("Number not found in the line.")
                break
            if l1t in line:
                found_target = True

        found_target = False
        lines_read = 0
        skip_counter = 0
        lines_to_read = 4 * (math.ceil(natoms / 6)) - 1
        start_index = 0
        for line in file:
            if found_target and lines_read < lines_to_read:
                if skip_counter != 3:  # Skip every fourth line
                    data_line = line.strip()
                    lines_read += 1
                    start_index += 1
                    parts = data_line.split("  ")
                    print(parts)
                    m = start_index
                    for j in range(1, len(parts)):
                        f[m - 1] = float(parts[j])
                        m = m + 3
                if skip_counter == 3:
                    lines_read += 1
                    start_index = 15 + start_index
                skip_counter = (skip_counter + 1) % 4
            if l2t in line:
                found_target = True
                print(line)
                # Skip the header line
                next(file)

        f = -f
        f = np.where(f == -0.0, 0.0, f)
        
        # Update the forces in the Molecule object
        molecule.update_forces(f)
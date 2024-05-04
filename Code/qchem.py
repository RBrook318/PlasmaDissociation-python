# Qchem.py 15/11/2023
import init
import subprocess
import os
import numpy as np 
import re
import math 
np.set_printoptions(precision =30)

def file_contains_string(file_path, search_string):
    with open(file_path, "r") as file:
        for line in file:
            if search_string in line:
                return True
    return False

def create_qchem_input(output_file, molecule, scf_algorithm="DIIS", Guess=True):
    # Extract molecule information
    symbols = molecule.symbols
    coordinates_bohr = molecule.coordinates  # Coordinates are in Bohr
    dissociation_flags = molecule.dissociation_flags
    multiplicity = molecule.multiplicity
    # Convert coordinates from Bohr to Angstrom
    bohr_to_angstrom = 0.529177
    coordinates_angstrom = coordinates_bohr * bohr_to_angstrom

    # Filter indices based on dissociation flag
    active_indices = [i for i, flag in enumerate(dissociation_flags) if flag == 'NO']
    # Q-Chem input file content
    qchem_input = (
        "$molecule\n"
        f"0 {multiplicity}\n"
        + "".join([f"{symbols[i]}   {coordinates_angstrom[i, 0]:.12f}   {coordinates_angstrom[i, 1]:.12f}   {coordinates_angstrom[i, 2]:.12f}\n" for i in active_indices])
        + "$end\n"
        "$rem\n"
        "    JOBTYPE             Force\n"
        "    EXCHANGE            BHHLYP\n"
        "    BASIS               6-31+G*\n"
        "    UNRESTRICTED        True\n"
        "    MAX_SCF_CYCLES      500\n"
        "    SYM_IGNORE          True\n"
        f"    SCF_Algorithm       {scf_algorithm}\n"  # Use the specified SCF algorithm
        "\n"
        "    SPIN_FLIP           True\n"
        "    SET_Iter            500\n"
        "\n"
        "    MAX_CIS_CYCLES      500\n"
    )

    # Add SCF_GUESS line if Guess is True
    if Guess:
        qchem_input += "    SCF_GUESS           Read\n"

    qchem_input += (
        "\n"
        "CIS_N_ROOTS 1\n"
        "CIS_STATE_DERIV 1\n"
        "$end\n"
    )
    # Write the Q-Chem input content to the output file
    with open(output_file, "w") as qchem_file:
        qchem_file.write(qchem_input)

def run_qchem(ncpu,file, molecule, n, nstates, Guess=True): 

    def submit_qchem_job(ncpu,file,output):
        subprocess.run(["qchem", "-save", "-nt", str(ncpu), file, output, "wf"])

    # Prepare f.inp file
    create_qchem_input(file, molecule, scf_algorithm="DIIS", Guess=Guess)

    # Submit the initial QChem job
    submit_qchem_job(ncpu,file,"f.out")



    if file_contains_string("f.out", "Thank you very much for using Q-Chem"):
        # Job completed successfully
        readqchem('f.out', molecule, n, nstates)
        pass   
    else:
        # Retry with a different setup
        create_qchem_input("f.inp", molecule, scf_algorithm="DIIS_GDM", Guess=False)
        submit_qchem_job(ncpu, "f.inp", "f2.out")

        # Check for the second failure
        if file_contains_string("f2.out", "Thank you very much for using Q-chem"):
            # Job failed both times, print error and exit
            readqchem('f2.out', molecule, n, nstates)
        else:
            with open("ERROR", "w") as file:
                file.write("Error occurred during QChem job.\n" + os.getcwd())
           

    # Append f.out content to f.all (is this super necessary??)
    with open("f.out", "r") as f_out, open("f.all", "a") as f_all:
        f_all.write(f_out.read())

def readqchem(output_file, molecule, natoms, nst):

    reduced_natoms = sum(flag.lower() != 'yes' for flag in molecule.dissociation_flags)
    ndim = 3 * reduced_natoms
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
                    print("energy = ", e)
                else:
                    print("Number not found in the line.")
                break
            if l1t in line:
                found_target = True

        found_target = False
        lines_read = 0
        skip_counter = 0
        lines_to_read = 4 * (math.ceil(natoms / 6)) - 1
        print("lines = ", lines_to_read)
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
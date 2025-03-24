# Qchem.py 15/11/2023
"""
QChem.py

This module provides functions to set up, run, and parse electronic structure calculations using QChem. It includes 
support for standard and spin-flip SCF calculations, geometry optimization, and normal mode analysis. This file is 
primarily used in conjunction with molecular dynamics simulations that require quantum chemical force calculations.

Dependencies:
- pyqchem
- numpy
- math
- os

Functions:
----------
file_contains_string(file_path: str, search_string: str) -> bool
    Checks if a specified string exists within a given file.

create_qchem_input(molecule: Structure, spin_flip: int, scf_algorithm: str = "DIIS", Guess: bool = True) -> QchemInput
    Generates a QChem input object based on molecule data, SCF algorithm, and spin-flip specification.

run_qchem(ncpu: int, molecule: Structure, n: int, nstates: int, spin_flip: int, Guess: bool = True) -> Structure
    Executes QChem calculation, updates molecule with electronic structure coefficients, and handles potential SCF 
    algorithm errors by retrying with alternative settings.

readqchem(output: str, molecule: Structure, natoms: int, nst: int, spin_flip: int) -> None
    Parses QChem output for total energy and forces, updating the molecule's energy and force attributes.

initial_conditions(symbols: list, coords: np.ndarray, cores: int) -> tuple
    Performs geometry optimization and normal mode analysis, outputting optimized geometry and frequency data.

Usage:
------
This module is intended to be imported into larger simulation scripts where QChem-based calculations are needed.
To perform an electronic structure calculation, set up a `molecule` object with required properties and call
`run_qchem` with the desired options.


"""
import os
import numpy as np 
import math
from pyqchem import QchemInput, Structure
from pyqchem import get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_frequencies import basic_frequencies
from pyqchem.parsers.basic import basic_parser_qchem
import subprocess
import time
import global_vars as gv

np.set_printoptions(precision =30)

# Electronic structure via QChem.
        
def file_contains_string(file_path, search_string):
    """
    Checks if a specified string exists within a file.

    Parameters
    ----------
    file_path : str
        Path to the file to be checked.
    search_string : str
        The string to search for within the file.

    Returns
    -------
    bool
        True if the search_string is found in the file, False otherwise.
    """
    with open(file_path, "r") as file:
        for line in file:
            if search_string in line:
                return True
    return False

def create_qchem_nac_input(molecule, scf_algorithm="DIIS", Guess=True):
    """
    Generates a QChem input configuration for simultaneous force and NAC calculations.

    Parameters
    ----------
    molecule : Structure
        A molecule structure with coordinates, symbols, and multiplicity set.
    scf_algorithm : str, optional
        The SCF convergence algorithm, by default "DIIS".
    Guess : bool, optional
        Whether to include initial SCF coefficients as the guess, by default True.


    Returns
    -------
    QchemInput
        A QChem input object with settings for simultaneous force and NAC calculations.
    """

    # Filter indices based on dissociation flag (if applicable)
    active_indices = [i for i, flag in enumerate(molecule.dissociation_flags) if flag == 'NO']
    active_coords = [molecule.coordinates[i] for i in active_indices]
    active_symbols = [molecule.symbols[i] for i in active_indices]


    # Create the molecule structure
    molecule = Structure(coordinates=active_coords, symbols=active_symbols, multiplicity=molecule.multiplicity)
    
    qc_inp = QchemInput(
        molecule,
        exchange='BHHLYP',
        basis=gv.basis,
        unrestricted=True,
        max_scf_cycles=500,
        sym_ignore=True,
        scf_algorithm=scf_algorithm,
        cis_n_roots=gv.num_states,
        extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500,'calc_nac':'true','CIS_DER_NUMSTATE': total_states},
        )

    
    # Get the Q-Chem input file content
    qc_content = qc_inp.get_txt().splitlines()  

    # Modify the content if Guess is True
    if Guess:
        for idx, line in enumerate(qc_content):
            if "scf_algorithm" in line:
                qc_content.insert(idx + 1, "SCF_GUESS           Read")
                break  # Stop after inserting once

    # Write the modified content back
    with open("ec.inp", "w") as f:
        f.write("\n".join(qc_content) + "\n")  # Ensure newlines are preserved
        f.write("\n$derivative_coupling\n")
        f.write("comment\n")
        for i in range(1, gv.num_states + 1):
            for j in range(i + 1, gv.num_states + 1):
                f.write(f"{i} {j}\n")
        f.write("$end\n")

    return qc_inp

def create_qchemforces_input(molecule, scf_algorithm="DIIS", Guess=True, excited_state=1):
    """
    Generates a QChem input configuration for an SCF calculation based on 
    molecular data, SCF algorithm, and spin-flip configuration.

    Parameters
    ----------
    molecule : Structure
        A molecule structure with coordinates, symbols, and multiplicity set.
    spin_flip : int
        Indicates whether to use spin-flip (1) or standard (0) DFT calculations.
    scf_algorithm : str, optional
        The SCF convergence algorithm, by default "DIIS".
    Guess : bool, optional
        Whether to include initial SCF coefficients as the guess, by default True.
    excited_state : int, optional
        The excited state for which gradients will be calculated, by default 1 (first excited state).
    basis : str, optional
        The basis set to use, by default '6-31+G*'.
    exchange : str, optional
        The exchange-correlation functional, by default 'BHHLYP'.

    Returns
    -------
    QchemInput
        A QChem input object with settings for the requested calculation type.

    Notes
    -----
    - When spin_flip is set to 1, additional parameters for spin-flip are 
      included in the input. If `Guess` is True, initial SCF coefficients 
      are included in the input for a better starting guess.
    - The `excited_state` parameter defines which excited state is used in 
      derivative calculations, affecting `cis_state_deriv`.
    """
    # Validate inputs
    if gv.spin_flip not in [0, 1]:
        raise ValueError("spin_flip must be 0 (standard) or 1 (spin-flip).")
    if excited_state < 1:
        raise ValueError("excited_state must be at least 1.")

    # Filter indices based on dissociation flag
    active_indices = [i for i, flag in enumerate(molecule.dissociation_flags) if flag == 'NO']
    active_coords = [molecule.coordinates[i] for i in active_indices]
    active_symbols = [molecule.symbols[i] for i in active_indices]
    coefficients = molecule.elecinfo

    # Create the molecule structure
    molecule = Structure(coordinates=active_coords, symbols=active_symbols, multiplicity=molecule.multiplicity)

    # Define common parameters
    if gv.spin_flip == 0:
        if Guess:
            qc_inp = QchemInput(
                molecule,
                jobtype='force',
                exchange='BHHLYP',
                scf_guess = coefficients,
                basis=gv.basis,
                unrestricted=True,
                max_scf_cycles=500,
                sym_ignore=True,
                scf_algorithm=scf_algorithm,
                extra_rem_keywords={'input_bohr': 'true'}
            )
        else:
            qc_inp = QchemInput(
                molecule,
                jobtype='force',
                exchange='BHHLYP',
                basis=gv.basis,
                unrestricted=True,
                max_scf_cycles=500,
                sym_ignore=True,
                scf_algorithm=scf_algorithm,
                extra_rem_keywords={'input_bohr': 'true'}
            )
    elif gv.spin_flip == 1: 
        if Guess:
            qc_inp = QchemInput(
                molecule,
                jobtype='force',
                exchange='BHHLYP',
                scf_guess = coefficients,
                basis=gv.basis,
                unrestricted=True,
                max_scf_cycles=500,
                sym_ignore=True,
                scf_algorithm=scf_algorithm,
                cis_n_roots=gv.num_states,
                cis_state_deriv= excited_state,
                extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500}
            )
        else:
            qc_inp = QchemInput(
                molecule,
                jobtype='force',
                exchange='BHHLYP',
                basis=gv.basis,
                unrestricted=True,
                max_scf_cycles=500,
                sym_ignore=True,
                scf_algorithm=scf_algorithm,
                cis_n_roots=gv.num_states,
                cis_state_deriv= excited_state,
                extra_rem_keywords={'input_bohr':'true','spin_flip':'true','set_iter':500}
            )

    # Add SCF guess if requested
    # if Guess:
    #     qc_inp.update_input({
    #         'scf_guess': coefficients
    #     })

    return qc_inp 
            
def run_qchem(molecule, Guess=True): 
    """
    Executes a QChem job to calculate the electronic structure for a given molecule
    and updates its properties based on the computed results.

    Parameters
    ----------
    ncpu : int
        The number of processors to use for the QChem calculation.
    molecule : Structure
        The molecular structure object with attributes such as coordinates, symbols, 
        and dissociation flags.
    Guess : bool, optional
        Whether to include an initial guess for SCF coefficients, by default True.

    Returns
    -------
    Structure
        The molecule object with updated electronic properties and forces.

    Notes
    -----
    - The function first attempts a QChem calculation using the 'DIIS' algorithm.
      If this fails, it retries with the 'DIIS_GDM' algorithm.
    - If both attempts fail, an error message is logged to the "ERROR" file, and the 
      program exits.
    - Upon successful completion, forces and energy data are extracted from the QChem 
      output, which is then passed to the `readqchem` function for detailed parsing.
    """
    if gv.num_states == 1:
        store_full_output = False
    else:
        store_full_output = True
        
    for state in range(1,gv.num_states+1):
        qc_inp=create_qchemforces_input(molecule,Guess=Guess,excited_state=state)
        try:
            qtime1 = time.time()
            output, ee = get_output_from_qchem(qc_inp,processors=gv.cores,return_electronic_structure=True,store_full_output = store_full_output)
            molecule.elecinfo=(ee['coefficients'])
            qtime2= time.time()
            print("time for qchem inside:",qtime2-qtime1)
            molecule.time[0] += qtime2-qtime1
            molecule.time[2] += qtime2-qtime1
            # Append f.out content to f.all
            with open("f.all", "a") as f_all:
                f_all.write(output)
        except:
            # print('Using DIIS_GDM algorithm')
            # Retry with a different setup
            
            qc_inp=create_qchemforces_input(molecule, scf_algorithm="DIIS_GDM", Guess=False,excited_state=state)
            try:
                qtime2= time.time()
                molecule.time[0] += qtime2-qtime1
                molecule.time[2] += qtime2-qtime1
                qtime1 = time.time()
                output, ee = get_output_from_qchem(qc_inp,processors=gv.cores,return_electronic_structure=True,store_full_output=store_full_output)
                molecule.elecinfo=(ee['coefficients'])
                qtime2= time.time()
                molecule.time[0] += qtime2-qtime1
                molecule.time[2] += qtime2-qtime1
                with open("f.all", "a") as f_all:
                    f_all.write(output)
            except:
                with open("ERROR", "w") as file:
                    file.write("Error occurred during QChem job. Help.\n" + os.getcwd())
                    file.write(output)
                exit()
        # Job completed successfully
        readqchemforces(output, molecule,state)

    if gv.num_states > 1:
        # First attempt using DIIS
        create_qchem_nac_input(molecule, scf_algorithm="DIIS", Guess=Guess)
        qtime1 = time.time()
        subprocess.run(["qchem", "-save", "-nt", str(gv.cores), "ec.inp", "ec.out", "wf"])
        qtime2 = time.time()
        molecule.time[0] += qtime2 - qtime1
        molecule.time[2] += qtime2 - qtime1

        try:
            readqchemnac('ec.out', molecule, gv.num_states)
        except ValueError:
            # If NAC section not found, retry with DIIS_GDM
            print("Derivative coupling not found, retrying with DIIS_GDM...")
            create_qchem_nac_input(molecule, scf_algorithm="DIIS_GDM", Guess=False)
            qtime1 = time.time()
            subprocess.run(["qchem", "-save", "-nt", str(gv.cores), "ec.inp", "ec.out", "wf"])
            qtime2 = time.time()
            molecule.time[0] += qtime2 - qtime1
            molecule.time[2] += qtime2 - qtime1
            
            try:
                readqchemnac('ec.out', molecule)
            except ValueError:
                raise ValueError("SF-CIS derivative coupling not found even after retrying with DIIS_GDM.")

        with open("ec.out", "r") as ec_file, open("ec.all", "a") as ec_all:
            ec_all.write(ec_file.read())
            ec_all.write("---------------------------------------------------\n")


    return molecule

def readqchemforces(output, molecule, state):

    """
    Parses QChem output to extract the SCF energy and gradient forces for a given molecule,
    and updates the molecule’s attributes accordingly.

    Parameters
    ----------
    output : str
        The raw output string from a QChem calculation.
    molecule : Structure
        The molecule object whose SCF energy and forces are to be updated.
    nst : int
        The number of states being considered.


    Notes
    -----
    - Extracts either the SCF energy or state-specific energy from the QChem output 
      depending on the calculation type.
    - Parses the gradient forces and updates the molecule's forces attribute.
    - Forces are formatted to avoid any negative zeros by replacing them with zeros.
    """
    reduced_natoms = sum(flag.lower() != 'yes' for flag in molecule.dissociation_flags)
    ndim = 3 * reduced_natoms
    if gv.spin_flip==1:
        enum = output.find(f'Total energy for state  {state}:')
        scf_erg = float(output[enum: enum+100].split()[5])
        molecule.scf_energy[state-1]=scf_erg
        l2t = ' Gradient of the state energy (including CIS Excitation Energy)'
    else:
        enum = output.find('Total energy in the final basis set')
        scf_erg = float(output[enum: enum+100].split()[8])
        molecule.scf_energy[state-1]=scf_erg
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
    f = f.reshape(-1,3)
    
    # Update the forces in the Molecule object
    # molecule.update_forces(f)
    molecule.forces[:,:,state-1] = f 
 
def readqchemnac(filename, molecule):
    """
    Parses QChem output to extract the SF-CIS derivative coupling for a given molecule,
    ensuring the data is structured correctly as (natoms, 3).

    Parameters
    ----------
    filename : str
        Path to the QChem output file.
    molecule : Structure
        The molecule object whose derivative coupling is to be updated.


    Returns
    -------
    np.ndarray
        An array of shape (natoms, 3) containing the derivative coupling values.

    Notes
    -----
    - Extracts the SF-CIS derivative coupling data from the QChem output.
    - Ensures the extracted data is correctly reshaped to (natoms, 3).
    """
    natoms = len(molecule.symbols)  # Number of atoms
    coupling_data = []
    found_coupling = False

    with open(filename, 'r') as f:
        for line in f:
            if "SF-CIS derivative coupling with ETF" in line:
                found_coupling = True
                next(f)  # Skip the header line
                next(f)
                continue
            
            if found_coupling:
                if line.strip() == "":
                    break  # Stop at an empty line
                parts = line.split()
                if len(parts) == 4:  # Expecting format: atom_index x y z
                    _, x, y, z = parts
                    coupling_data.append([float(x), float(y), float(z)])

    if not found_coupling:
        raise ValueError("SF-CIS derivative coupling section not found in the output.")

    # Convert to numpy array
    coupling_array = np.array(coupling_data, dtype=float)

    # Ensure correct shape (natoms, 3)
    if coupling_array.shape != (natoms, 3):
        raise ValueError(f"Incorrect derivative coupling shape: {coupling_array.shape}, expected ({natoms}, 3)")
    molecule.coupling = coupling_array

    return coupling_array

def initial_conditions(symbols,coords,cores):
    """
    Generates initial conditions for a molecule by optimizing its geometry 
    and calculating vibrational modes using QChem.

    Parameters
    ----------
    symbols : list of str
        Atomic symbols of the atoms in the molecule.
    coords : list of list of float
        Atomic coordinates for each atom in the molecule.
    cores : int
        The number of cores to use for QChem calculations.

    Returns
    -------
    opt_geoms : Structure
        The optimized molecular geometry as a Structure object.
    parser_output : dict
        Parsed vibrational frequency data from the QChem frequency calculation.

    Notes
    -----
    - First runs a geometry optimization calculation for the molecule, saving output to 
      'optimisation.out'.
    - Then performs a vibrational frequency calculation and saves the output to 'modes.out'.
    - Utilizes basic_optimization and basic_frequencies parsers to interpret QChem output.
    """
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
    opt_coords=(parser_output['optimized_molecule'])
    values=opt_coords.get_coordinates()
    atom=opt_coords.get_symbols()
    opt_geoms = ''
    for i in range(len(atom)):
        opt_geoms+=(atom[i]+" "+str(float(values[i][0]))+" "+str(float(values[i][1]))+" "+str(float(values[i][2]))+" \n")
    # Find Normal modes
    qc_inp = QchemInput(opt_coords,
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



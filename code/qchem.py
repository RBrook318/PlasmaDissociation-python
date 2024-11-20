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

def create_qchem_input(molecule, spin_flip,scf_algorithm="DIIS", Guess=True):
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

    Returns
    -------
    QchemInput
        A QChem input object with settings for the requested calculation type.

    Notes
    -----
    - When spin_flip is set to 1, additional parameters for spin-flip are 
      included in the input. If `Guess` is True, initial SCF coefficients 
      are included in the input for a better starting guess.
    """
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
                      
def run_qchem(ncpu, molecule, nstates, spin_flip, Guess=True): 
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
    nstates : int
        Number of electronic states to consider in the calculation.
    spin_flip : int
        Indicates whether spin-flip (1) or standard (0) calculations are used.
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
    readqchem(output, molecule,nstates,spin_flip)
    # Append f.out content to f.all
    # with open("f.all", "a") as f_all:
    #     f_all.write(output)
    return molecule


def readqchem(output, molecule, nst,spin_flip):

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
    spin_flip : int
        Indicates if spin-flip calculations (1) or standard (0) calculations are used.

    Notes
    -----
    - Extracts either the SCF energy or state-specific energy from the QChem output 
      depending on the calculation type.
    - Parses the gradient forces and updates the molecule's forces attribute.
    - Forces are formatted to avoid any negative zeros by replacing them with zeros.
    """
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



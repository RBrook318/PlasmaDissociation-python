# This should be the code that handles the multiple different electronic structure methods that can be used
# Current valid methods are QChem and pySCF(?)
import os
import numpy as np 
import re
import math
from pyqchem import QchemInput, Structure
from pyqchem import get_output_from_qchem
from pyscf import gto, scf, dft, grad

np.set_printoptions(precision =30)


def run_elec_structure(molecule, ncpu,n,nstates,spin_flip, method,Guess):
    if method == 'QChem':
        run_qchem(ncpu, molecule,n, nstates,spin_flip,Guess=Guess)
    elif method == 'PySCF': 
        run_pySCF(molecule,ncpu)

# -------------------------------------------------------------------
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
       qc_inp.update_input({'scf_guess':coefficients})  
 
    return qc_inp
                      
def run_qchem(ncpu, molecule, n, nstates, spin_flip, Guess=True): 
    qc_inp=create_qchem_input(molecule, spin_flip, scf_algorithm="DIIS", Guess=Guess)
    try:
        output, ee = get_output_from_qchem(qc_inp,processors=ncpu,return_electronic_structure=True)
        molecule.update_elecinfo(ee['coefficients'])
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
    
# -------------------------------------------------------------------------
# PySCF section 

import init 
from pyscf import gto, scf, dft, grad

def create_pyscf_molecule(molecule):
    # Extract molecule information
    symbols = molecule.symbols
    coordinates_bohr = molecule.coordinates  # Coordinates are in Bohr
    dissociation_flags = molecule.dissociation_flags
    multiplicity = molecule.multiplicity
    
    # Filter indices based on dissociation flag
    active_indices = [i for i, flag in enumerate(dissociation_flags) if flag == 'NO']
    
    mol = gto.Mole()
    mol.atom = [(symbols[i], coordinates_bohr[i]) for i in active_indices]
    mol.basis = '6-31+g*'
    mol.spin = multiplicity - 1
    mol.charge = 0
    mol.unit = 'Bohr'  # Set unit to Bohr
    mol.build()

    return mol

def run_pyscf_calculation(mol, scf_algorithm, exchange_functional=None):
    if exchange_functional:
        if mol.spin == 0:
            mf = dft.RKS(mol)
        else:
            mf = dft.UKS(mol)
        # Manually specify BHHLYP components
        if exchange_functional.upper() == 'BHHLYP':
            mf.xc = '0.5*HF + 0.5*B88,LYP'
        else:
            mf.xc = exchange_functional
    else:
        if scf_algorithm == "DIIS":
            mf = scf.RHF(mol) if mol.spin == 0 else scf.UHF(mol)
        elif scf_algorithm == "DIIS_GDM":
            mf = scf.newton(mol)
        else:
            raise ValueError(f"Unknown SCF algorithm: {scf_algorithm}")

    mf.conv_tol = 1e-10
    mf.max_cycle = 500
    
    energy = mf.kernel()
    
    # Calculate forces (gradients)
    if mol.spin == 0:
        if exchange_functional:
            grad_calc = grad.RKS(mf)
        else:
            grad_calc = grad.RHF(mf)
    else:
        if exchange_functional:
            grad_calc = grad.UKS(mf)
        else:
            grad_calc = grad.UHF(mf)
    
    forces = grad_calc.kernel()
    
    return energy, forces

def run_pySCF(molecule, npcu, exchange_functional='BHHLYP'):
    e = molecule.scf_energy

    mol = create_pyscf_molecule(molecule)
    energy, forces = run_pyscf_calculation(mol, 'DIIS', exchange_functional)
    e[1] = energy
    molecule.update_scf_energy(e)
    molecule.update_forces(-forces)
    print(molecule.forces)


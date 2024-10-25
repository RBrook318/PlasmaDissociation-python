# PySCF section 
import time
import os
from pyscf import gto, scf, dft, grad, lib, hessian
# from pyscf.geomopt.geometric_solver import optimize
# from pyscf.prop.freq import rks as freq
from gpu4pyscf.dft import rks, uks

def create_pyscf_molecule(molecule):
    """
    Create a PySCF molecule object based on the provided molecule data.
    :param molecule: A molecule object containing symbols, coordinates in Bohr, dissociation flags, etc.
    :return: A PySCF molecule object.
    """
    symbols = molecule.symbols
    coordinates_bohr = molecule.coordinates  # Coordinates in Bohr
    dissociation_flags = molecule.dissociation_flags
    multiplicity = molecule.multiplicity
    
    # Only include atoms with 'NO' dissociation flags
    active_indices = [i for i, flag in enumerate(dissociation_flags) if flag == 'NO']
    
    mol = gto.Mole()
    mol.atom = [(symbols[i], coordinates_bohr[i]) for i in active_indices]
    mol.basis = '6-31+G*'  # Pre-set basis
    mol.spin = multiplicity-1
    mol.charge = 0
    mol.unit = 'Bohr'
    mol.build()

    return mol

def run_pyscf_calculation(mol, scf_algorithm='DIIS', prev_mf=None):
    """
    Run a PySCF SCF calculation on the given molecule.
    :param mol: A PySCF molecule object.
    :param scf_algorithm: SCF algorithm to use (default: 'DIIS').
    :param prev_mf: A previous mean-field solution to start from, if available.
    :return: A tuple containing (SCF energy, forces, mean-field object).
    """

    # Set the appropriate SCF method based on the spin
    if mol.spin == 0:
        mf = dft.RKS(mol).density_fit()
    else:
        mf = dft.UKS(mol).density_fit()

    mf.xc = '0.5*HF + 0.5*B88,LYP'
    mf.conv_tol = 1e-8
    mf.conv_tol_grad = 1e-7
    mf.max_cycle = 500
    mf.verbose = 4
    
    # Restart from a previous mean-field solution if available
    if prev_mf is not None:
        mf.mo_coeff = prev_mf.mo_coeff
        mf.mo_occ = prev_mf.mo_occ

    # Run SCF energy calculation and measure time
    timescf1 = time.time()
    energy = mf.kernel()
    timescf2 = time.time()
    print('SCF energy time:', timescf2 - timescf1)
    
    # Compute gradients and measure time
    timescf1 = time.time()
    if mol.spin == 0:
        grad_calc = grad.RKS(mf)
    else:
        grad_calc = grad.UKS(mf)
    forces = grad_calc.kernel()
    timescf2 = time.time()
    print('SCF gradients time:', timescf2 - timescf1)

    return energy, forces, mf
 

def run_pySCF(molecule,Guess=True,use_gpu = False):
    e = molecule.scf_energy
    mol = create_pyscf_molecule(molecule)
    guess = False
    time1 = time.time()
    if use_gpu==False:
        if Guess ==False: 
            energy, forces,mf = run_pyscf_calculation(mol, 'DIIS')
        else:
            energy, forces,mf = run_pyscf_calculation(mol, 'DIIS',molecule.elecinfo)
    else: 
        if Guess ==False: 
            energy, forces,mf = run_pyscf_calculation_gpu(mol, 'DIIS')
        else:
            energy, forces,mf = run_pyscf_calculation_gpu(mol, 'DIIS',molecule.elecinfo)        
    time2=time.time()
    print('pyscf job: ',  time2-time1)
    e[0] = energy
    molecule.update_scf_energy(e)
    forces = -forces
    forces = forces.reshape(-1) 
    molecule.update_forces(forces)
    molecule.update_elecinfo(mf)

def initial_conditions(symbols, coords):
    # Create PySCF molecule
    mol = gto.Mole()
    mol.atom = [(symbols[i], coords[i]) for i in range(len(symbols))]
    mol.basis = '6-31+G*'
    mol.unit = 'Bohr'
    mol.build()
    
    # Optimize geometry
    start_time = time.time()
    mf = dft.RKS(mol)
    mf.xc = '0.5*HF + 0.5*B88,LYP'
    mf.conv_tol = 1e-7
    mf.max_cycle = 500
    mf.kernel()
    mol_eq = optimize(mf)
    print(mol_eq.atom_coords())
    geom_opt = mol_eq.atom_coords()
    end_time = time.time()
    print(f'Geometry optimization took {end_time - start_time} seconds')
    
    # Calculate normal modes (frequencies)
    start_time = time.time()
    mf = dft.freq(mol_eq).run()
    w, modes = freq.Freq(mf).kernel()
    print(freq)
    end_time = time.time()
    print(f'Frequency calculation took {end_time - start_time} seconds')
    
    return geom_opt, freq

def run_pyscf_calculation_gpu(mol, scf_algorithm='DIIS', prev_mf=None):
    """
    Run a PySCF SCF calculation on the GPU.
    :param mol: A PySCF molecule object.
    :param scf_algorithm: SCF algorithm to use (default: 'DIIS').
    :param prev_mf: A previous mean-field solution to start from, if available.
    :return: A tuple containing (SCF energy, forces, mean-field object).
    """

    # Set the appropriate SCF method based on the spin and move to GPU
    if mol.spin == 0:
        mf = rks.RKS(mol).to_gpu()
    else:
        mf = uks.UKS(mol).to_gpu().density_fit() 
    mf.xc = '0.5*HF + 0.5*B88,LYP'
    mf.scf_algorithm = 'DIIS'
    mf.conv_tol = 1e-9
    mf.conv_tol_grad = 1e-9
    mf.max_cycle = 500
    mf.verbose = 3
    
    # Restart from a previous mean-field solution if available
    if prev_mf is not None:
        mf.mo_coeff = prev_mf.mo_coeff
        mf.mo_occ = prev_mf.mo_occ

    # Run SCF energy calculation and measure GPU time
    timescf1 = time.time()
    energy = mf.kernel()
    timescf2 = time.time()
    print('GPU SCF energy time:', timescf2 - timescf1)
    
    # Compute gradients and measure GPU time
    timescf1 = time.time()
    if mol.spin == 0:
        grad_calc = mf.nuc_grad_method().to_gpu()
    else:
        grad_calc = mf.nuc_grad_method().to_gpu()  # For UKS
    forces = grad_calc.kernel()
    timescf2 = time.time()
    print('GPU SCF gradients time:', timescf2 - timescf1)

    return energy, forces, mf
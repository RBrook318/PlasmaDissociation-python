# PySCF section 
import time
import os
from pyscf import gto, scf, dft, grad, lib, hessian
# from pyscf.geomopt.geometric_solver import optimize
# from pyscf.prop.freq import rks as freq
# from gpu4pyscf.dft import rks

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

def run_pyscf_calculation(mol, scf_algorithm, prev_mf=None,use_gpu =False):
    
    print(lib.num_threads())

    if mol.spin == 0:
        mf = dft.RKS(mol).density_fit()
    else:
        mf = dft.UKS(mol).density_fit()
    
    mf.xc = '0.5*HF + 0.5*B88,LYP'
    mf.conv_tol = 1e-8
    mf.conv_tol_grad = 1e-7
    mf.max_cycle = 100
    mf.verbose = 3

    if prev_mf is not None:
        mf.mo_coeff = prev_mf.mo_coeff
        mf.mo_occ = prev_mf.mo_occ
    if use_gpu:
        mf = rks.RKS(mol).to_gpu()
        mf.xc = '0.5*HF + 0.5*B88,LYP'
        mf = rks.RKS(mol).run()
        gobj = mf.nuc_grad_method().to_gpu()
        forces = gobj.kernel()
        energy = mf.to_gpu().kernel()  
    else: 
        timescf1= time.time()
        energy = mf.kernel()
        timescf2= time.time()
        print('SCF energy time:', timescf2-timescf1)
        timescf1= time.time()
        if mol.spin == 0:
            grad_calc = grad.RKS(mf)
        else:
            grad_calc = grad.UKS(mf)
        forces = grad_calc.kernel()
        timescf2= time.time()
        print('SCF gradients time:', timescf2-timescf1)
    return energy, forces, mf 

def run_pySCF(molecule,Guess=True,use_gpu = False):
    e = molecule.scf_energy
    mol = create_pyscf_molecule(molecule)
    guess = False
    time1 = time.time()
    if Guess ==False: 
        energy, forces,mf = run_pyscf_calculation(mol, 'DIIS')
    else:
        energy, forces,mf = run_pyscf_calculation(mol, 'DIIS',molecule.elecinfo)
    time2=time.time()
    print('pyscf job: ',  time2-time1)
    e[1] = energy
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

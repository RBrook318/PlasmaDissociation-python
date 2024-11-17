# PySCF section 
import time
import os
from pyscf import gto, scf, dft, grad, lib, hessian
from pyscf.geomopt.geometric_solver import optimize
# from pyscf.prop.freq import rks as freq


def create_pyscf_molecule(molecule):
    """
    Creates a PySCF molecule object based on the provided molecule data.

    Parameters
    ----------
    molecule : Structure
        A molecule object containing symbols, coordinates in Bohr, dissociation flags, and multiplicity.

    Returns
    -------
    mol : PySCF gto.Mole
        A PySCF molecule object with active atoms included based on dissociation flags.

    Notes
    -----
    - Only includes atoms with dissociation flags set to 'NO' for active calculations.
    - The PySCF molecule is built with '6-31+G*' basis and default Bohr units.
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
    Runs a PySCF SCF calculation on the given molecule.

    Parameters
    ----------
    mol : gto.Mole
        A PySCF molecule object.
    scf_algorithm : str, optional
        SCF algorithm to use, default is 'DIIS'.
    prev_mf : SCF object, optional
        Previous mean-field solution for initialization if available.

    Returns
    -------
    energy : float
        The calculated SCF energy.
    forces : ndarray
        The calculated nuclear gradients (forces).
    mf : SCF object
        The mean-field SCF object with solution information.

    Notes
    -----
    - Uses RKS for closed-shell systems and UKS for open-shell.
    - Supports density fitting and hybrid functional setup for computational efficiency.
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
    mf.verbose = 3
    
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
    """
    Runs a PySCF calculation with optional GPU support and updates the molecule object with results.

    Parameters
    ----------
    molecule : Structure
        The molecule object to update with SCF energy, forces, and electronic structure information.
    Guess : bool, optional
        If True, uses a previous guess for SCF initialization.
    use_gpu : bool, optional
        If True, runs the calculation using GPU support.

    Notes
    -----
    - Uses `run_pyscf_calculation` or `run_pyscf_calculation_gpu` based on the GPU flag.
    - Updates molecule's SCF energy, forces, and electronic information attributes.
    """
    e = molecule.scf_energy
    mol = create_pyscf_molecule(molecule)

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
    """
    Initializes a molecule with optimized geometry and frequency data using PySCF.

    Parameters
    ----------
    symbols : list of str
        Atomic symbols of the atoms in the molecule.
    coords : list of list of float
        Atomic coordinates in Bohr units.

    Returns
    -------
    geom_opt : ndarray
        The optimized atomic coordinates.
    frequencies : ndarray
        Vibrational frequencies (in atomic units).
    normal_modes : ndarray
        Normal mode displacement vectors.

    Notes
    -----
    - First, performs a geometry optimization and then calculates vibrational frequencies.
    - Utilizes RKS method with hybrid functional '0.5*HF + 0.5*B88,LYP' for both steps.
    """
    # Create PySCF molecule
    mol = gto.Mole()
    mol.atom = [(symbols[i], coords[i]) for i in range(len(symbols))]
    mol.basis = '6-31+G*'
    mol.unit = 'Bohr'
    mol.build()

    # Optimize geometry
    start_time = time.time()
    mf = dft.RKS(mol)
    mf.verbose = 4
    mf.xc = '0.5*HF + 0.5*B88,LYP'
    mf.conv_tol = 1e-7
    mf.max_cycle = 500
    mf.kernel()
    mol_eq = optimize(mf)
    geom_opt = mol_eq.atom_coords()
    end_time = time.time()
    print(f'Geometry optimization took {end_time - start_time} seconds')

    # Calculate normal modes (frequencies)
    start_time = time.time()
    freq_calc = freq.Freq(mol_eq, mf)
    frequencies, normal_modes = freq_calc.kernel()
    end_time = time.time()
    print(f'Frequency calculation took {end_time - start_time} seconds')

    # Output results
    print("Optimized Geometry:")
    print(geom_opt)
    print("Frequencies (a.u.):")
    print(frequencies)
    print("Normal Modes:")
    print(normal_modes)

    return geom_opt, frequencies, normal_modes

def run_pyscf_calculation_gpu(mol, scf_algorithm='DIIS', prev_mf=None):
    from gpu4pyscf.dft import rks, uks
    """
    Runs a GPU-accelerated PySCF SCF calculation.

    Parameters
    ----------
    mol : gto.Mole
        A PySCF molecule object.
    scf_algorithm : str, optional
        SCF algorithm to use, default is 'DIIS'.
    prev_mf : SCF object, optional
        Previous mean-field solution for initialization if available.

    Returns
    -------
    energy : float
        The calculated SCF energy.
    forces : ndarray
        The calculated nuclear gradients (forces).
    mf : SCF object
        The mean-field SCF object with solution information.

    Notes
    -----
    - Uses GPU-accelerated RKS/UKS for closed/open-shell calculations.
    - Supports density fitting with hybrid functionals for enhanced performance.
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
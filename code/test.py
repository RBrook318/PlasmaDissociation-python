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
    mol.basis = '6-31+g*' #def2-TZVP
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
    molecule.update_forces(forces)

molecule = init.initialize_structure(2,0)

run_pySCF(molecule,8)
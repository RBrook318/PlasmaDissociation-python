# This should be the code that handles the multiple different electronic structure methods that can be used
# Current valid methods are QChem and pySCF(?)




# Currently supported electronic structure methods: QChem, PySCF, and GPUPySCF

def run_elec_structure(molecule, ncpu,nstates, spin_flip, method, Guess,basis):
    """
    Executes an electronic structure calculation on a molecule using a specified method.

    Parameters
    ----------
    molecule : Molecule
        The molecular structure and properties.
    ncpu : int
        Number of CPUs to use for the calculation.
    n : int
        Calculation identifier or iteration.
    nstates : int
        Number of electronic states to consider.
    spin_flip : bool
        Toggle spin-flip calculations.
    method : str
        The electronic structure method to use ('QChem', 'PySCF', or 'GPUPySCF').
    Guess : object
        Initial guess for electronic structure calculations.

    Returns
    -------
    molecule : Molecule
        Updated molecule object with new electronic structure data.

    Raises
    ------
    ValueError
        If an unsupported method is specified.
    """
    if method == 'QChem':
        import qchem as qc
        qc.run_qchem(ncpu, molecule, nstates, spin_flip,basis=basis, Guess=Guess,)
    elif method == 'PySCF': 
        import pyscf as pyscf
        pyscf.run_pySCF(molecule, Guess)
    elif method == 'GPUPySCF':
        import pyscf as pyscf
        pyscf.run_pySCF(molecule, Guess, use_gpu=True)
    else:
        raise ValueError(f"Unsupported method: {method}. Choose from 'QChem', 'PySCF', or 'GPUPySCF'.")

    return molecule
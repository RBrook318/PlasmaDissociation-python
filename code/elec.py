# This should be the code that handles the multiple different electronic structure methods that can be used
# Current valid methods are QChem and pySCF(?)




def run_elec_structure(molecule, ncpu,n,nstates,spin_flip, method,Guess):
    if method == 'QChem':
        import qchem as qc
        qc.run_qchem(ncpu, molecule,n, nstates,spin_flip,Guess=Guess)
    elif method == 'PySCF': 
        import elec_pyscf as pyscf
        pyscf.run_pySCF(molecule,Guess)
    elif method == 'GPUPySCF':
        import elec_pyscf as pyscf
        pyscf.run_pySCF(molecule,Guess,use_gpu=True)
    return molecule
        
# -------------------------------------------------------------------

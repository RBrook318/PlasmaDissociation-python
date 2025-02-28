"""
Propagator Module for Molecular Dynamics Simulations

This module contains functions for propagating the molecular system over time,
including applying Magnus expansion, calculating force vectors, and handling dissociation events.

Functions:
---------
- CompForceEhr: Computes the Ehrenfest force vector based on amplitudes, forces, energies, and couplings.
- magnus2: Computes the Magnus expansion for a given Hamiltonian and time step.
- exp_matrix: Evaluates the matrix exponential, supporting a scaling factor.
- calculate_electronic_hamiltonian: Computes the electronic Hamiltonian matrix based on molecule energy and velocities.
- shrink_molecule: Reduces molecule size by removing dissociated atoms.
- restore_molecule: Restores the molecule to its original form after processing.
- prop_1: Propagates the molecular system using a modified first-order Magnus expansion method.
- prop_2: Propagates the molecular system with a second-order Magnus expansion, adjusting energy and forces iteratively.
- fragments: Identifies fragments of the molecule and adjusts the spin if necessary.
- prop_diss: Propagates dissociated atoms by moving them based on their velocities.

"""

import numpy as np
from scipy.linalg import expm
import init


np.set_printoptions(precision=30)

def CompForceEhr(A, F, E, C, nst):
    """
    Computes the Ehrenfest force vector using state amplitudes, forces, and couplings.

    Parameters:
    - A : (nst,) ndarray
        Array of amplitudes.
    - F : (natoms, 3, nst) ndarray
        Array of forces for each atom in each state.
    - E : (nst,) ndarray
        Array of energies for each state.
    - C : (natoms, 3) ndarray
        Coupling values for each atom.
    - nst : int
        Number of states.

    Returns:
    - ForceVector : (natoms, 3) ndarray
        Calculated Ehrenfest force vector.
    """
    natoms, _, _ = F.shape
    ForceVector = np.zeros((natoms, 3), dtype=np.float64)
    f1 = np.zeros((natoms, 3), dtype=np.float64)
    f2 = np.zeros((natoms, 3), dtype=np.float64)

    # f1 calculation (weighted force sum)
    if nst == 1:
        f1 += F[:, :, 0] * np.abs(A[0])**2
    else:
        for i in range(nst):
            f1 += F[:, :, i] * np.abs(A[i])**2

    # f2 calculation (non-adiabatic coupling term)
    if nst > 1:
        for i in range(nst):
            for j in range(i+1, nst):
                ae = 2.0 * np.real(np.conj(A[i]) * A[j]) * (E[i] - E[j])
                f2 += ae * C[:,:,i,j]

    ForceVector = f1 + f2
    
    return ForceVector

def magnus2(H0, H1, dt):
    """
    Calculate the Magnus expansion for the given Hamiltonians H0 and H1 at a given time step dt.

    Parameters:
    - H0, H1: 2D NumPy arrays representing the Hamiltonians.
    - dt: Time step.

    Functions:
    - exp_matrix - prop.py

    Returns:
    - magH: Resulting Magnus expansion matrix.
    """
    ndim = H0.shape[0]

    # Calculate the average
    # print("electronic hamilonians \n")
    # print(H0,"\n")
    # print(H1,"\n")
    # print("Average: \n")

    Hav = np.sum(np.diag(H0) + np.diag(H1)) / (2 * ndim)
    # print(Hav, "\n")
    # print(Hav)
    # The trace matrix
    Htr = np.diag(np.full(ndim, Hav))
    # print("Trace matrix \n")
    # print(Htr)
    
    a0 = (H1 + H0) / 2.0 - Htr
    W1 = dt * a0
    
    # Assuming exp_pade is a function that performs the matrix exponential using Pade approximation
    magH = exp_matrix(W1) * np.exp(Hav * dt)
    # print("MagH with timestep",dt,":")
    # print(magH)
    # print("----------")
    return magH

def exp_matrix(A, t_in=None):
    """
    Computes the matrix exponential, optionally scaled by a given factor.

    Parameters:
    ----------
    - A : np.ndarray
        The square matrix to exponentiate.
    - t_in : float, optional
        Scaling factor for the exponentiation. Defaults to 1.0 if not provided.

    Returns:
    -------
    - expA : np.ndarray
        The matrix exponential of A scaled by t_in.

    Notes:
    -----
    Uses the SciPy `expm` function unless the input matrix is near zero, 
    in which case an identity matrix is returned.
    """
    m = A.shape[0]
    expA = np.zeros((m, m), dtype=np.complex256)

    t = 1.0
    if t_in is not None:
        t = t_in

    if np.max(np.abs(A)) < 1.0e-7:
        expA = np.eye(m, dtype=np.complex256)
    else:
        expA = expm(A * t)


    return expA

def calculate_electronic_hamiltonian(molecule, velocities, coupling):
    """
    Calculates the electronic Hamiltonian for a molecular system, 
    incorporating energy shifts and velocity-dependent couplings.

    Parameters:
    ----------
    - molecule : object
        Molecular object containing SCF energy values.
    - velocities : np.ndarray
        Array of velocities for each atom.
    - coupling : np.ndarray
        Coupling matrix between states.

    Returns:
    -------
    - electronic_hamiltonian : np.ndarray
        The computed electronic Hamiltonian matrix.

    Notes:
    -----
    The diagonal elements are the SCF energies shifted by a constant (77.67785291), 
    while the off-diagonal elements represent velocity-coupling terms.
    """
    # print("Calculating the hamiltonian")
    nst = len(molecule.scf_energy)
    natoms = len(molecule.symbols)
    # print("Number of states:", nst)
    ii = 1j  # imaginary unit

    electronic_hamiltonian = np.zeros((nst, nst), dtype=np.cdouble)

    for n1 in range(nst):
        electronic_hamiltonian[n1, n1] = molecule.scf_energy[n1] + 77.67785291 
        for n2 in range(n1 + 1, nst):
            electronic_hamiltonian[n1, n2] = -ii * np.sum(velocities * coupling[:,:,n1,n2])
            electronic_hamiltonian[n2, n1] = -electronic_hamiltonian[n1, n2]

    return electronic_hamiltonian

def shrink_molecule(molecule):
    """
    Reduces a molecule by removing atoms marked for dissociation, 
    returning a copy of the molecule with only the remaining atoms.

    Parameters:
    ----------
    - molecule : object (molecule class)
        Molecular object containing attributes such as symbols, coordinates, and momenta,
        with dissociation flags marking atoms to be removed.

    Returns:
    -------
    - shrunk_molecule : object (molecule class)
        Copy of the original molecule with dissociated atoms removed.
    - shrunk_index : list of int
        List of indices corresponding to retained atoms.

    Notes:
    -----
    This function creates a copy of the input molecule and removes atoms flagged 
    for dissociation. The `shrunk_index` can later be used to map back to the 
    original molecule's indexing.
    """
    natoms = len(molecule.symbols)
    dis_index = []
    shrunk_index = []

    for i in range(natoms): 
        if molecule.dissociation_flags[i] == 'YES':
            dis_index.append(i)
        else:
            shrunk_index.append(i)

    shrunk_molecule = molecule.copy()
    
    new_symbols = [molecule.symbols[i] for i in range(natoms) if i not in dis_index]
    new_coordinates = [molecule.coordinates[i,:] for i in range(natoms) if i not in dis_index]
    new_momenta = [molecule.momenta[i,:] for i in range(natoms) if i not in dis_index]

    # Update the new molecule with the modified information
    shrunk_molecule.update_symbols(new_symbols)
    shrunk_molecule.update_coordinates(new_coordinates)
    shrunk_molecule.update_momenta(new_momenta)
    
    
    return shrunk_molecule, shrunk_index

def restore_molecule(molecule, shrunk_molecule, shrunk_index):

    """
    Restores atoms and properties to a molecule from a shrunk version, 
    based on the specified index mapping. Used at the end of a prop step 
    to update a full molecule after it was shrunk.

    Parameters:
    ----------
    - molecule : object
        The original molecule object to be restored, containing all atoms and properties.
    - shrunk_molecule : object
        The modified molecule object with dissociated atoms removed.
    - shrunk_index : list of int
        List of indices mapping the shrunk molecule's atoms back to the original molecule.

    Returns:
    -------
    - molecule : object
        The restored molecule with updated symbols, coordinates, momenta, and other properties.

    Raises:
    ------
    - ValueError
        If any index in `shrunk_index` is out of bounds for the original molecule.

    Notes:
    -----
    This function updates the molecule's symbols, coordinates, and momenta using 
    data from the shrunk molecule based on `shrunk_index`. Other properties such 
    as SCF energy, forces, amplitudes, and timestep are also restored.
    """
    for i, index in enumerate(shrunk_index):
        if index < 0 or index >= len(molecule.symbols):
            raise ValueError(f"Invalid index: {index}")
        # Update attributes based on shrunk index

        molecule.symbols[index] = shrunk_molecule.symbols[i]
        molecule.coordinates[index,:] = shrunk_molecule.coordinates[i,:]
        molecule.momenta[index,:] = shrunk_molecule.momenta[i,:]
        
    molecule.update_scf_energy(shrunk_molecule.scf_energy)
    molecule.update_forces(shrunk_molecule.forces)
    molecule.update_amplitudes(shrunk_molecule.amplitudes)
    molecule.update_timestep(shrunk_molecule.timestep)

    return molecule

def prop_1(molecule1, molecule2, natoms, nst, increment):
    """
    Propagates molecular amplitudes, coordinates, and momenta over a time increment,
    considering forces and velocities based on current state information.

    Parameters:
    ----------
    - molecule1 : object
        Input molecule object representing the initial molecular state.
    - molecule2 : object
        Second molecule object to receive updated properties after propagation.
    - natoms : int
        Number of atoms in the molecule.
    - nst : int
        Number of states in the electronic Hamiltonian.
    - increment : float
        Time step over which propagation occurs.

    Returns:
    -------
    - molecule2 : object
        Updated molecule object with propagated properties.

    Functions:
    ---------
    - `shrink_molecule(molecule)` : -prop.py
        Shrinks a molecule to exclude atoms marked for dissociation.

    - `restore_molecule(molecule, shrunk_molecule, shrunk_index)` : -prop.py
        Restores the molecule to its original full form after propagation.

    - `calculate_electronic_hamiltonian(shrunk_molecule, velocities, Coupling)` : - prop.py
        Calculates the electronic Hamiltonian.

    - `magnus2(ham_1, ham_2, increment)` : - prop.py
        Calculates time evolution of electronic amplitudes using the Magnus expansion.

    - `CompForceEhr(amplitudes, forces, scf_energy, Coupling, nst)` : -prop.py
        Computes the Ehrenfest forces for the system.

    Notes:
    -----
    This function calculates molecular velocities, electronic Hamiltonian, and 
    force vectors to update the amplitudes and forces over a set time increment. 
    Uses intermediate computations for accurate propagation, storing results in 
    temporary molecule `molecule2`.
    """
    amplitudes = molecule1.amplitudes
    velocities = np.zeros((natoms,3))
    forces_1 = molecule1.forces
    scf_energy_1 = molecule1.scf_energy
    molecule2 = molecule1.copy()

    shrunk_molecule, shrunk_index = shrink_molecule(molecule1)
    
    natoms = len(shrunk_molecule.symbols)
    
    coupling = np.zeros((natoms,3,nst,nst))
    for i in range(0, natoms):
        velocities[i,:] = shrunk_molecule.momenta[i,:]/shrunk_molecule.masses[i]
        if nst > 1:
            for y in range(0,3):
                coupling[i,y,0,1] = shrunk_molecule.coupling[i,y]  
                coupling[i,y,1,0] = -shrunk_molecule.coupling[i,y]  
            
    
    Eham_1 = calculate_electronic_hamiltonian(shrunk_molecule, velocities, coupling)
    # print("magnus2 output shape:", magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20).shape)
    print(Eham_1)
    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment/20), amplitudes)
    # print('Amplitudes:', Amplitudes_temp)
    Force_vector=CompForceEhr(amplitudes,forces_1,scf_energy_1,coupling,nst)/10
    
    for im in range(1, 10):
        
        A1 = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment/10), Amplitudes_temp)
        # print('Amplitudes:', A1)
        Amplitudes_temp = A1
        Force_vector +=  CompForceEhr(Amplitudes_temp, forces_1, scf_energy_1, coupling,nst)/10
    
    print("Final amplitudes of prop 1: \n")
    print(A1)
    print("------------------\n")

    shrunk_molecule.update_amplitudes(A1)
    shrunk_molecule.update_timestep(shrunk_molecule.timestep+increment)
  
    for i in range(natoms):
        shrunk_molecule.coordinates[i,:] = shrunk_molecule.coordinates[i,:] + increment*velocities[i,:] + ((increment**2)/2)*Force_vector[i,:]/shrunk_molecule.masses[i]
        shrunk_molecule.momenta[i,:] = shrunk_molecule.momenta[i,:] + increment*Force_vector[i,:]
    restore_molecule(molecule2, shrunk_molecule, shrunk_index)

    return molecule2

def prop_2(molecule1, molecule2, natoms, nst, increment):
    """
    Interpolates properties between two molecular states and propagates 
    amplitudes, energy, and forces based on electronic Hamiltonians and couplings.

    Parameters:
    ----------
    - molecule1 : object
        Initial molecular state to propagate.
    - molecule2 : object
        Secondary molecule object used for interpolation.
    - natoms : int
        Number of atoms in the molecule.
    - nst : int
        Number of states in the electronic Hamiltonian.
    - increment : float
        Time step over which propagation occurs.

    Returns:
    -------
    - molecule1 : object
        Updated molecule object with propagated properties.

    Functions:
    ---------
    - `shrink_molecule(molecule)` : - prop.py
        Shrinks a molecule to exclude atoms marked for dissociation.

    - `restore_molecule(molecule, shrunk_molecule, shrunk_index)` : - prop.py
        Restores the molecule to its original full form after propagation.

    - `calculate_electronic_hamiltonian(shrunk_molecule, velocities, Coupling)` : -prop.py
        Calculates the electronic Hamiltonian.

    - `magnus2(ham_1, ham_2, increment)` : - prop.py
        Calculates time evolution of electronic amplitudes using the Magnus expansion.

    - `CompForceEhr(amplitudes, forces, scf_energy, Coupling, nst)` : - prop.py
        Computes the Ehrenfest forces for the system.

    Notes:
    -----
    Uses the shrunk versions of the molecules to calculate interpolated velocities, 
    energy, forces, and electronic Hamiltonians. The force vector is averaged over 
    intermediate steps and applied to propagate the state of `molecule1`.
    """




    shrunk_molecule1, shrunk_index = shrink_molecule(molecule1)
    shrunk_molecule2, shrunk_index = shrink_molecule(molecule2)



    natoms = len(shrunk_molecule1.symbols)
    velocities_1 = np.zeros((natoms,3))
    velocities_2 = np.zeros((natoms,3))
    coupling_1 = np.zeros((natoms,3,nst,nst))
    coupling_2 = np.zeros((natoms,3,nst,nst))
    for i in range(0, natoms):
        velocities_1[i,:] = shrunk_molecule1.momenta[i,:]/shrunk_molecule1.masses[i]
        velocities_2[i,:] = shrunk_molecule2.momenta[i,:]/shrunk_molecule2.masses[i]
        if nst >1: 
            for y in range(0,3):
                coupling_1[i,y,0,1] = shrunk_molecule1.coupling[i,y]  
                coupling_1[i,y,1,0] = -shrunk_molecule1.coupling[i,y]  
                coupling_2[i,y,0,1] = shrunk_molecule2.coupling[i,y]  
                coupling_2[i,y,1,0] = -shrunk_molecule2.coupling[i,y]  

    Eham_1 = calculate_electronic_hamiltonian(shrunk_molecule1,velocities_1,coupling_1)
    # print("Coupling difference:", coupling_1 - coupling_2)
    Eham_2 = calculate_electronic_hamiltonian(shrunk_molecule2,velocities_2,coupling_2)
    # print(Eham_1)
    # print(Eham_2)

    ElPhase = np.ones(nst)  # Initialize ElPhase[1] = 1 (Fortran index 1 → Python index 0)
    for j in range(1, nst):  
        num = np.sum(coupling_1[:, :, 0, j] * coupling_2[:, :, 0, j])  # Selecting (natoms,) elements
        den = np.sqrt(np.sum(coupling_1[:, :, 0, j]**2) * np.sum(coupling_2[:, :, 0, j]**2))
    
        val = num / den if den != 0 else 0.0  # Avoid division by zero
        ElPhase[j] = np.sign(val)  # Equivalent to sign(1., val)
        # print("Phase:", ElPhase[j])
        # Warning condition
        if abs(val) < 0.5 and abs(shrunk_molecule2.amplitudes[j]) >= 0.35:
            print(f"!! Warning: the sign for state {j} is not reliable! {val:.4f}")
    
    for i in range(1, nst):
        coupling_2[:,:,i,:] *= ElPhase[i]  # Apply phase to the second dimension
        coupling_2[:,:,:,i] *= ElPhase[i]  # Apply phase to the third dimension
    
    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20), shrunk_molecule1.amplitudes)
    Energy_temp = 0.05 * shrunk_molecule2.scf_energy + 0.95 * shrunk_molecule1.scf_energy
    Forces_temp = 0.05 * shrunk_molecule2.forces + 0.95 * shrunk_molecule1.forces
    Coupling_temp = 0.05 * coupling_2 + 0.95 * coupling_1

    Force_vector = CompForceEhr(Amplitudes_temp, Forces_temp, Energy_temp, Coupling_temp, nst)/10.
    for im in range(1, 10):
        Eham_temp = (im * Eham_2 + (10 - im) * Eham_1) * 0.1
        Energy_temp = ((0.1 * im) + 0.05) * shrunk_molecule2.scf_energy + (0.95 - im * 0.1) * shrunk_molecule1.scf_energy
        Forces_temp = ((0.1 * im) + 0.05) * shrunk_molecule2.forces + (0.95 - im * 0.1) * shrunk_molecule1.forces
        Coupling_temp = ((0.1 * im) + 0.05) * coupling_2 + (0.95 - im * 0.1) * coupling_1
        A1 = np.matmul(magnus2(-1j * Eham_temp, -1j * Eham_temp, increment/10), Amplitudes_temp)
        Amplitudes_temp = A1
        # print("Calculated amplitudes: \n")
        # print(A1)
        # print("------------------\n")
        Force_vector_temp = CompForceEhr(A1, Forces_temp, Energy_temp, Coupling_temp, nst)
        Force_vector += Force_vector_temp / 10


    A1 = np.matmul(magnus2(-1j * Eham_2, -1j * Eham_2, increment / 20), A1)
    print("Final amplitudes: \n")
    print(A1)
    print("------------------\n")
    # print(shrunk_molecule1.momenta)
    # print('Force_vector: ',Force_vector)
    shrunk_molecule1.momenta = shrunk_molecule1.momenta + increment * Force_vector 
    # print(shrunk_molecule1.momenta)
    
    shrunk_molecule1.update_symbols(shrunk_molecule2.symbols)
    shrunk_molecule1.update_coordinates(shrunk_molecule2.coordinates)
    shrunk_molecule1.update_scf_energy(shrunk_molecule2.scf_energy)
    shrunk_molecule1.update_forces(shrunk_molecule2.forces)
    shrunk_molecule1.update_amplitudes(A1)
    shrunk_molecule1.update_timestep(shrunk_molecule2.timestep)
    restore_molecule(molecule1,shrunk_molecule1, shrunk_index)
    molecule1.coupling = molecule2.coupling

    return molecule1

def fragments(molecule, spin_flip):
    """
    Analyzes and marks atoms in a molecule that meet dissociation criteria and
    updates multiplicity based on the presence of a spin flip.

    Parameters:
    ----------
    - molecule : object
        The molecular object representing atomic and molecular states.
    - spin_flip : int
        Indicates whether spin flip is to be applied: 
        1 to increase multiplicity, 0 to decrease if certain conditions are met.

    Returns:
    -------
    - molecule : object
        Updated molecule object with dissociation flags set for atoms meeting
        dissociation criteria.
    - dissociated : int
        Indicator of whether any atoms were marked as dissociated (1 for true, 0 for false).

    Notes:
    -----
    - Iterates over atoms with no dissociation flags, checking if their 
      forces are below a defined threshold.
    - Updates the molecular multiplicity and applies a spin flip if required.
    - Now operates with forces stored in a `(natoms, 3, nstates)` array.
    """
    natom = len(molecule.symbols)
    molecules_with_no_flag = [i for i in range(natom) if molecule.dissociation_flags[i] == 'NO']

    dissociated = 0
    j = 1
    to_remove = []  # Collect indices to remove

    for i in molecules_with_no_flag:
        # Compute force magnitude across all states
        val = np.sqrt(np.sum(molecule.forces[j-1, :, :]**2))
        if val < 1.0e-5:
            molecule.multiplicity -= 1
            if spin_flip == 1 and molecule.multiplicity == 2:
                molecule.multiplicity = 4
            elif spin_flip == 0 and molecule.multiplicity == 0:
                molecule.multiplicity = 2

            molecule.dissociation_flags[i] = 'YES'
            dissociated = 1
            to_remove.append(j-1)
        j += 1
    # Remove dissociated atoms from forces
    if to_remove:
        molecule.forces = np.delete(molecule.forces, to_remove, axis=0)
        molecule.coupling = np.delete(molecule.coupling, to_remove,axis=0)

    return molecule, dissociated

def prop_diss(molecule, increment): 
    """
    Propagates coordinates of dissociated atoms based on their velocities, updating their
    positions over the specified time increment.

    Parameters:
    ----------
    - molecule : object
        The molecular object, containing current atom positions, velocities, and masses.
    - increment : float
        Time step used to propagate dissociated atoms’ coordinates.

    Returns:
    -------
    - molecule : object
        Updated molecule object with propagated coordinates for dissociated atoms.


    Notes:
    -----
    This function calculates the velocity of atoms marked as dissociated and uses these
    to propagate their positions over the time increment. Positions are updated only for
    atoms with dissociation flags set to 'YES'.
    """
    natoms = len(molecule.symbols)
    dis_index = []
    shrunk_index = []


    for i in range(natoms): 
        if molecule.dissociation_flags[i] == 'YES':
            dis_index.append(i)
        else:
            shrunk_index.append(i)

    
    velocities = np.zeros((3))

    for i in dis_index:
        velocities[:] = molecule.momenta[i, :] / molecule.masses[i]
        molecule.coordinates[i, :] = molecule.coordinates[i, :] + increment * velocities[:]

    return molecule
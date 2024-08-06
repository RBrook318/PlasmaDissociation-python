import numpy as np
from scipy.linalg import expm


np.set_printoptions(precision=30)

def CompForceEhr(A, F, E, C,nst):
    ndim = F.shape
    ForceVector = np.zeros(ndim, dtype=np.float64)
    f1 = np.zeros(ndim, dtype=np.float64)
    f2 = np.zeros(ndim, dtype=np.float64)

    for i in range(nst-1):
        f1 += F[:] * np.abs(A[i])**2
    
    
    for i in range(nst):
        for j in range(i+1, nst):
            ae = 2.0 * np.real(np.conj(A[i]) * A[j]) * (E[i] - E[j])
            f2 += ae * C
   
    ForceVector = f1 + f2

    return ForceVector

def magnus2(H0, H1, dt):
    """
    Calculate the Magnus expansion for the given Hamiltonians H0 and H1 at a given time step dt.

    Parameters:
    - H0, H1: 2D NumPy arrays representing the Hamiltonians.
    - dt: Time step.

    Returns:
    - magH: Resulting Magnus expansion matrix.
    """
    ndim = H0.shape[0]

    # Calculate the average
    Hav = np.sum(np.diag(H0) + np.diag(H1)) / (2 * ndim)
    
    # The trace matrix
    Htr = np.diag(np.full(ndim, Hav))

    
    a0 = (H1 + H0) / 2.0 - Htr
    W1 = dt * a0
    
    # Assuming exp_pade is a function that performs the matrix exponential using Pade approximation
    magH = exp_matrix(W1) * np.exp(Hav * dt)

    return magH

def exp_matrix(A, t_in=None):
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
    nst = len(molecule.scf_energy)
    ii = 1j  # imaginary unit

    electronic_hamiltonian = np.zeros((nst, nst), dtype=np.cdouble)

    for n1 in range(nst):
        electronic_hamiltonian[n1, n1] = molecule.scf_energy[n1] + 77.67785291 
        for n2 in range(n1 + 1, nst):
            electronic_hamiltonian[n1, n2] = -ii * np.sum(velocities * coupling)
            electronic_hamiltonian[n2, n1] = -electronic_hamiltonian[n1, n2]

    return electronic_hamiltonian

def shrink_molecule(molecule):

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
    # Make sure the lengths match

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
    amplitudes = molecule1.amplitudes
    velocities = np.zeros((natoms,3))
    forces_1 = molecule1.forces
    scf_energy_1 = molecule1.scf_energy
    molecule2 = molecule1.copy()

    shrunk_molecule, shrunk_index = shrink_molecule(molecule1)
    
    natoms = len(shrunk_molecule.symbols)
    
    Coupling = 0
    for i in range(0, natoms):
        velocities[i,:] = shrunk_molecule.momenta[i,:]/shrunk_molecule.masses[i]
    
    Eham_1 = calculate_electronic_hamiltonian(shrunk_molecule, velocities, Coupling)

    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20), amplitudes.reshape(-1, 1))

    Force_vector=CompForceEhr(amplitudes,forces_1,scf_energy_1,Coupling,nst)/10
    
    for im in range(1, 10):
        A1 = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 10), Amplitudes_temp)
        Amplitudes_temp = A1
       
        Force_vector +=  CompForceEhr(Amplitudes_temp, forces_1, scf_energy_1, Coupling,nst)/10
        

    shrunk_molecule.update_amplitudes(A1)
    shrunk_molecule.update_timestep(shrunk_molecule.timestep+increment)
  
    Force_vector = Force_vector.reshape(-1, 3) 
    for i in range(natoms):
        shrunk_molecule.coordinates[i,:] = shrunk_molecule.coordinates[i,:] + increment*velocities[i,:] + ((increment**2)/2)*Force_vector[i,:]/shrunk_molecule.masses[i]
        shrunk_molecule.momenta[i,:] = shrunk_molecule.momenta[i,:] + increment*Force_vector[i,:]

    restore_molecule(molecule2, shrunk_molecule, shrunk_index)

    return molecule2

def prop_2(molecule1, molecule2, natoms, nst, increment):
    ndim = natoms*3
    velocities_1 = np.zeros((natoms,3))
    velocities_2 = np.zeros((natoms,3))



    shrunk_molecule1, shrunk_index = shrink_molecule(molecule1)
    shrunk_molecule2, shrunk_index = shrink_molecule(molecule2)


    Coupling = 0
    natoms = len(shrunk_molecule1.symbols)
    for i in range(0, natoms):
        velocities_1[i,:] = shrunk_molecule1.momenta[i,:]/shrunk_molecule1.masses[i]
        velocities_2[i,:] = shrunk_molecule2.momenta[i,:]/shrunk_molecule2.masses[i]
        
    Eham_1 = calculate_electronic_hamiltonian(shrunk_molecule1,velocities_1,Coupling)

    Eham_2 = calculate_electronic_hamiltonian(shrunk_molecule2,velocities_2,Coupling)

    
    ElPhase = np.zeros(nst, dtype=int)
    ElPhase[0] = 1

    for j in range(1, nst):
        # val = np.sum(Coupling * Coupling) / np.sqrt(np.sum(Coupling**2) * np.sum(Coupling**2))
        val = 0
        ElPhase[j] = np.where(val >= 0, 1, -1)
        if abs(val) < 0.5 and abs(shrunk_molecule2.amplitudes[j]) >= 0.35:
            print(f'!! Warning: the sign for state {j} is not reliable! {val:.4f}')

    for i in range(1, nst):
        Coupling *= ElPhase[i]
        Coupling *= ElPhase[i]

    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20), shrunk_molecule1.amplitudes)
    Energy_temp = 0.05 * shrunk_molecule2.scf_energy + 0.95 * shrunk_molecule1.scf_energy
    Forces_temp = 0.05 * shrunk_molecule2.forces + 0.95 * shrunk_molecule1.forces
    Coupling_temp = 0.05 * Coupling + 0.95 * Coupling

    Force_vector = CompForceEhr(Amplitudes_temp, Forces_temp, Energy_temp, Coupling_temp, nst)/10.
    for im in range(1, 10):
        Eham_temp = (im * Eham_2 + (10 - im) * Eham_1) * 0.1
        Energy_temp = (0.1 * im + 0.05) * shrunk_molecule2.scf_energy + (0.95 - im * 0.1) * shrunk_molecule1.scf_energy
        Forces_temp = (0.1 * im + 0.05) * shrunk_molecule2.forces + (0.95 - im * 0.1) * shrunk_molecule1.forces
        Coupling_temp = (0.1 * im + 0.05) * Coupling + (0.95 - im * 0.1) * Coupling
        A1 = np.matmul(magnus2(-1j * Eham_temp, -1j * Eham_temp, increment / 10), Amplitudes_temp)
        Amplitudes_temp = A1
        Force_vector_temp = CompForceEhr(A1, Forces_temp, Energy_temp, Coupling_temp, nst)
        Force_vector += Force_vector_temp / 10

    A1 = np.matmul(magnus2(-1j * Eham_2, -1j * Eham_2, increment / 20), Amplitudes_temp)
    Force_vector = Force_vector.reshape(-1, 3)
    
    shrunk_molecule1.momenta = shrunk_molecule1.momenta + increment * Force_vector 
    
    shrunk_molecule1.update_symbols(shrunk_molecule2.symbols)
    shrunk_molecule1.update_coordinates(shrunk_molecule2.coordinates)
    shrunk_molecule1.update_scf_energy(shrunk_molecule2.scf_energy)
    shrunk_molecule1.update_forces(shrunk_molecule2.forces)
    shrunk_molecule1.update_amplitudes(A1)
    shrunk_molecule1.update_timestep(shrunk_molecule2.timestep)
    restore_molecule(molecule1,shrunk_molecule1, shrunk_index)
    

    return molecule1

def fragements(molecule,spin_flip):
    natom = len(molecule.symbols)
    molecules_with_no_flag = [i for i in range(1, natom + 1) if molecule.dissociation_flags[i - 1] == 'NO']

    dissociated = 0    
    j = 1

    for i in molecules_with_no_flag:
        val = np.sqrt(np.sum(molecule.forces[3 * (j - 1):3 * j]**2))

        if val < 1.0e-5:
 
            molecule.multiplicity -= 1
            if spin_flip ==1:
                if molecule.multiplicity == 2:
                    molecule.multiplicity = 4
            elif spin_flip ==0:
                if molecule.multiplicity == 0:
                    molecule.multiplicity = 2
            molecule.dissociation_flags[i - 1] = 'YES'
            dissociated = 1
            molecule.forces = np.delete(molecule.forces, 3 * (j-1))
            molecule.forces = np.delete(molecule.forces, 3 * (j-1))
            molecule.forces = np.delete(molecule.forces, 3 * (j-1))
            j -= 1
        j += 1

    return molecule, dissociated

def prop_diss(molecule, increment): 
    natoms = len(molecule.symbols)
    dis_index = []
    shrunk_index = []
    Mau = 1822.887

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


    
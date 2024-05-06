import numpy as np
from scipy.linalg import expm
import init
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
        print(molecule.scf_energy[n1], "+ 77.67785291 = ", electronic_hamiltonian[n1,n1])

        for n2 in range(n1 + 1, nst):
            electronic_hamiltonian[n1, n2] = -ii * np.sum(velocities * coupling)
            electronic_hamiltonian[n2, n1] = -electronic_hamiltonian[n1, n2]

    return electronic_hamiltonian

def prop_1(molecule1, molecule2, natoms, nst, increment):
    Mau=1822.887
    amplitudes = molecule1.amplitudes
    velocities = np.zeros((natoms,3))
    mass = np.zeros((natoms))
    forces_1 = molecule1.forces
    scf_energy_1 = molecule1.scf_energy
    molecule2 = molecule1.copy()
    Coupling = 0
    for i in range(0, natoms):
        if molecule1.symbols[i] == 'C':
            mass[i] = 12 * Mau
        elif molecule1.symbols[i] == 'N':
            mass[i] = 14 * Mau
        elif molecule1.symbols[i] == 'H':
            mass[i] = Mau
        elif molecule1.symbols[i] == 'D':
            mass[i] = 2 * Mau
        elif molecule1.symbols[i] == 'F':
            mass[i] = 19 * Mau
        elif molecule1.symbols[i] == 'O':
            mass[i] = 16 * Mau
        else:
            print('Atom', molecule1.symbols[i], 'is not supported')
            raise ValueError('Unsupported atom type')
        velocities[i,:] = molecule1.momenta[i,:]/mass[i]

    Eham_1 = calculate_electronic_hamiltonian(molecule1, velocities, Coupling)

    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20), amplitudes.reshape(-1, 1))

    Force_vector=CompForceEhr(amplitudes,forces_1,scf_energy_1,Coupling,nst)/10
    
    for im in range(1, 10):
        A1 = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 10), Amplitudes_temp)
        Amplitudes_temp = A1
       
        Force_vector +=  CompForceEhr(Amplitudes_temp, forces_1, scf_energy_1, Coupling,nst)/10
        

    molecule2.update_amplitudes(A1)
    molecule2.update_timestep(molecule1.timestep+increment)
  
    Force_vector = Force_vector.reshape(-1, 3) 

    for i in range(natoms):
        molecule2.coordinates[i,:] = molecule2.coordinates[i,:] + increment*velocities[i,:] + ((increment**2)/2)*Force_vector[i,:]/mass[i]

        molecule2.momenta[i,:] = molecule1.momenta[i,:] + increment*Force_vector[i,:]

    
    return molecule2

def prop_2(molecule1, molecule2, natoms, nst, increment):

    Mau=1822.887
    ndim = natoms*3
    velocities_1 = np.zeros((natoms,3))
    velocities_2 = np.zeros((natoms,3))
    mass = np.zeros((natoms))



    Coupling = 0
    natoms = len(molecule1.symbols)
    for i in range(0, natoms):
        if molecule1.symbols[i] == 'C':
            mass[i] = 12 * Mau
        elif molecule1.symbols[i] == 'N':
            mass[i] = 14 * Mau
        elif molecule1.symbols[i] == 'H':
            mass[i] = Mau
        elif molecule1.symbols[i] == 'D':
            mass[i] = 2 * Mau
        elif molecule1.symbols[i] == 'F':
            mass[i] = 19 * Mau
        else:
            print('Atom', molecule1.symbols[i], 'is not supported')
            raise ValueError('Unsupported atom type')
        velocities_1[i,:] = molecule1.momenta[i,:]/mass[i]
        velocities_2[i,:] = molecule2.momenta[i,:]/mass[i]
        
    Eham_1 = calculate_electronic_hamiltonian(molecule1,velocities_1,Coupling)

    Eham_2 = calculate_electronic_hamiltonian(molecule2,velocities_2,Coupling)

    
    ElPhase = np.zeros(nst, dtype=int)
    ElPhase[0] = 1

    for j in range(1, nst):
        # val = np.sum(Coupling * Coupling) / np.sqrt(np.sum(Coupling**2) * np.sum(Coupling**2))
        val = 0
        ElPhase[j] = np.where(val >= 0, 1, -1)
        if abs(val) < 0.5 and abs(molecule2.amplitudes[j]) >= 0.35:
            print(f'!! Warning: the sign for state {j} is not reliable! {val:.4f}')

    for i in range(1, nst):
        Coupling *= ElPhase[i]
        Coupling *= ElPhase[i]

    Amplitudes_temp = np.matmul(magnus2(-1j * Eham_1, -1j * Eham_1, increment / 20), molecule1.amplitudes)
    Energy_temp = 0.05 * molecule2.scf_energy + 0.95 * molecule1.scf_energy
    Forces_temp = 0.05 * molecule2.forces + 0.95 * molecule1.forces
    Coupling_temp = 0.05 * Coupling + 0.95 * Coupling

    Force_vector = CompForceEhr(Amplitudes_temp, Forces_temp, Energy_temp, Coupling_temp, nst)/10.

    for im in range(1, 10):
        Eham_temp = (im * Eham_2 + (10 - im) * Eham_1) * 0.1
        Energy_temp = (0.1 * im + 0.05) * molecule2.scf_energy + (0.95 - im * 0.1) * molecule1.scf_energy
        Forces_temp = (0.1 * im + 0.05) * molecule2.forces + (0.95 - im * 0.1) * molecule1.forces
        Coupling_temp = (0.1 * im + 0.05) * Coupling + (0.95 - im * 0.1) * Coupling
        A1 = np.matmul(magnus2(-1j * Eham_temp, -1j * Eham_temp, increment / 10), Amplitudes_temp)
        Amplitudes_temp = A1
        Force_vector_temp = CompForceEhr(A1, Forces_temp, Energy_temp, Coupling_temp, nst)
        Force_vector += Force_vector_temp / 10.
    A1 = np.matmul(magnus2(-1j * Eham_2, -1j * Eham_2, increment / 20), Amplitudes_temp)
    Force_vector = Force_vector.reshape(-1, 3)
    
    molecule1.momenta = molecule1.momenta + increment * Force_vector 
    
    molecule1.update_symbols(molecule2.symbols)
    molecule1.update_coordinates(molecule2.coordinates)
    molecule1.update_scf_energy(molecule2.scf_energy)
    molecule1.update_forces(molecule2.forces)
    molecule1.update_amplitudes(A1)
    molecule1.update_timestep(molecule2.timestep)


    return molecule1

# def fragements(molecule):
#     natom = len(molecule.symbols)
#     molecules_with_no_flag = [i for i in range(1, natom + 1) if molecule.dissociation_flags[i - 1] == 'NO']

#     dissociated = 0    
#     j = 1

#     for i in molecules_with_no_flag:
#         val = np.sqrt(np.sum(molecule.forces[3 * (j - 1):3 * j]**2))

#         if val < 1.0e-5:
 
#             molecule.multiplicity -= 1
#             if molecule.multiplicity == 1:
#                 molecule.multiplicity = 1
#             molecule.dissociation_flags[i - 1] = 'YES'
#             dissociated = 1
#             molecule.forces = np.delete(molecule.forces, 3 * (j-1))
#             molecule.forces = np.delete(molecule.forces, 3 * (j-1))
#             molecule.forces = np.delete(molecule.forces, 3 * (j-1))
#             j -= 1
#         j += 1

#     return molecule, dissociated

def prop_diss(molecule, increment): 
    dis_index = []
    Mau = 1822.887

    velocities = np.zeros((3))

    for i in dis_index:
        if molecule.symbols[i] == 'C':
            mass = 12 * Mau
        elif molecule.symbols[i] == 'N':
            mass = 14 * Mau
        elif molecule.symbols[i] == 'H':
            mass = Mau
        elif molecule.symbols[i] == 'D':
            mass = 2 * Mau
        elif molecule.symbols[i] == 'F':
            mass = 19 * Mau

        velocities[:] = molecule.momenta[1, :] / mass
        molecule.coordinates[1, :] = molecule.coordinates[1, :] + increment * velocities[:]

    return molecule

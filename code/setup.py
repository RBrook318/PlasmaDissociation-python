import numpy as np
import json
import pubchempy as pcp
import init
from scipy.spatial.distance import pdist, squareform
import os

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:,0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
    T=T*0.0000031668
    # Initialize arrays for random
    Meff = np.zeros(nmod)
    rv = np.zeros((nmod, mom_num))
    for i in range(nmod):
        for j in range(n):
            Meff[i] = Meff[i]+np.sum(((Ax[j, i]**2) + (Ay[j, i]**2) + (Az[j, i]**2)) * m[j])
        rv[i, :] = rn[i, :] * np.sqrt(2 * T / Meff[i])
    # Calculate the velocity by applying it through the tranformation matrix of normal modes.
    Vx = np.dot(Ax, rv)
    Vy = np.dot(Ay, rv)
    Vz = np.dot(Az, rv)
    Px = np.zeros((n,mom_num))
    Py = np.zeros((n,mom_num))
    Pz = np.zeros((n,mom_num))
    for i in range(n):
        Px[i,:] = Vx[i,:]*m[i]
        Py[i,:] = Vy[i,:]*m[i]
        Pz[i,:] = Vz[i,:]*m[i]
    
    return Px, Py, Pz

def get_geometry_file():

    symbols = []
    coordinates = []
    with open("initial_guess.txt", "r") as file:
        for line in file:
            parts = line.split()
            if len(parts) == 4:  # Ensuring correct format
                symbols.append(parts[0])  # Atomic symbol
                x, y, z = map(float, parts[1:])
                coordinates.append([x * 1.88973, y * 1.88973, z * 1.88973])  # Convert to Bohr


    return symbols, coordinates

import numpy as np



def momenta_checks(Px, Py, Pz, symbols, temperature, output_file="momentum_checks.txt"):
    """
    Calculate the kinetic energy for each repeat, compare it to the expected value,
    and calculate the average magnitude of momentum for each atom across all repeats.
    Results are saved to a file.

    Parameters:
        Px, Py, Pz (np.ndarray): Momenta in x, y, z directions (shape: [natoms, repeats]).
        symbols (list): List of atomic symbols (e.g., ['H', 'C', 'O']).
        temperature (float): Temperature in Kelvin.
        output_file (str): File to save the results.
    """
    natoms = Px.shape[0]  # Number of atoms
    repeats = Px.shape[1]  # Number of repeats

    # Get the masses of the atoms from the ATOMIC_MASSES dictionary
    masses = init.setup_masses(symbols)

    # Initialize variables to store kinetic energy for each repeat
    kinetic_energies = np.zeros(repeats)

    # Loop through each repeat
    for j in range(repeats):
        file_kinetic_energy = 0.0  # Initialize kinetic energy for the current repeat
        for i in range(natoms):
            # Calculate the magnitude of the momentum vector
            magnitude = np.sqrt(Px[i, j]**2 + Py[i, j]**2 + Pz[i, j]**2)

            # Calculate the kinetic energy for this atom
            kinetic_energy = 0.5 * (magnitude**2 / masses[i])

            # Add the kinetic energy to the total for this repeat
            file_kinetic_energy += kinetic_energy

        # Store the total kinetic energy for this repeat
        kinetic_energies[j] = file_kinetic_energy

    # Calculate the average kinetic energy across all repeats
    average_kinetic_energy = np.mean(kinetic_energies)

    # Calculate the expected kinetic energy based on the equation
    # Expected KE = (3N - 6) * 0.5 * kT, where k = 3.16881e-6 (atomic units)
    expected_ke = (3 * natoms - 6) * 3.16881e-6 * temperature

    # Calculate the average magnitude of momentum for each atom across all repeats
    avg_momentum_magnitude = np.zeros(natoms)
    for i in range(natoms):
        total_momentum_magnitude = 0.0
        for j in range(repeats):
            total_momentum_magnitude += np.sqrt(Px[i, j]**2 + Py[i, j]**2 + Pz[i, j]**2)
        avg_momentum_magnitude[i] = total_momentum_magnitude / repeats

    # Write results to a file
    with open(output_file, "w") as file:
        file.write("Kinetic Energy Checks:\n")
        file.write(f"Expected Kinetic Energy: {expected_ke:.6f}\n")
        file.write(f"Average Kinetic Energy: {average_kinetic_energy:.6f}\n")
        file.write("Average Momentum Magnitude per Atom:\n")
        for i in range(natoms):
            file.write(f"Atom {i+1} ({symbols[i]}): {avg_momentum_magnitude[i]:.6f}\n")
        for j in range(repeats):
            file.write(f"Repeat {j+1}:\n")
            file.write(f"  Kinetic Energy: {kinetic_energies[j]:.6f}\n")
            if np.isclose(kinetic_energies[j], expected_ke, rtol=0.1*expected_ke):
                file.write("  Kinetic Energy check passed.\n")
            else:
                file.write("  Kinetic Energy check failed.\n")
        file.write("\n")






def generate_bondarr(symbols, coordinates):

    BOND_THRESHOLDS = {
        ("C", "C"): [(2.7401, "C=C"), (3.25, "C-C")],  # 1.45 Å and 1.70 Å * 1.88973 for C=C and C-C
        ("C", "F"): [(2.8346, "C-F")],  # 1.50 Å * 1.88973
        ("C", "H"): [(2.1732, "C-H")],  # 1.15 Å * 1.88973
        ("C", "O"): [(2.75, "C-O")],  # 1.43 Å * 1.88973
        ("O", "H"): [(2.0, "O-H")],  # 1.00 Å * 1.88973
        ("H", "O"): [(2.0, "O-H")],  # 1.00 Å * 1.88973
    }   

    natoms = len(symbols)

    # Compute pairwise distances
    distances = squareform(pdist(coordinates))  # Already in Bohr

    with open("../results/bondarr.txt", "w") as file:
        for i in range(natoms):
            for j in range(i + 1, natoms):
                elem1, elem2 = symbols[i], symbols[j]
                d = distances[i, j]  # Distance between atom i and j in Bohr
                
                # Sort elements alphabetically to match dictionary keys
                key = tuple(sorted((elem1, elem2)))
        
                # Check if bond exists in dictionary
                if key in BOND_THRESHOLDS:

                    # Loop through possible bond types for this pair
                    for max_dist, bond_type in BOND_THRESHOLDS[key]:
                        if d <= max_dist:
                            # Write bond type to file
                            file.write(f"{i+1}-{j+1}:{bond_type}\n")
                            
    
    print("bondarr.txt successfully written!")


def get_geometry_pubchem(molecule_name):
    search = pcp.get_compounds(molecule_name, 'name', record_type='3d')
    molecule = search[0]
    bonds = molecule.bonds
    atoms = molecule.atoms
    
    # Create a mapping from atom IDs to symbols
    atom_id_to_symbol = {atom.aid: atom.element for atom in atoms}
    symbols = [atom.element for atom in atoms]
    coordinates = [[atom.x*1.88973, atom.y*1.88973, atom.z*1.88973] for atom in atoms]
    with open("../results/bondarr.txt", "w") as file:
        for bond in bonds:
            atom1 = bond.aid1
            atom2 = bond.aid2
            order = bond.order
            
            # Determine the bond type based on the order
            if order == 1:
                bond_type = "-"
            elif order == 2:
                bond_type = "="
            elif order == 3:
                bond_type = "#"
            else:
                bond_type = "-"  # Default to single bond if unknown order
            
            # Get the symbols for the atoms
            atom1_symbol = atom_id_to_symbol.get(atom1, "Unknown")
            atom2_symbol = atom_id_to_symbol.get(atom2, "Unknown")
            
            # Write to file in the format 1-3:C-H
            file.write(f"{atom1}-{atom2}:{atom1_symbol}{bond_type}{atom2_symbol}\n")

    return symbols,coordinates

def organise_modes(modes,atoms):
    numeric_modes = np.zeros((len(modes)*atoms, 3))
    cnt=0
    for i in range(len(modes)):
        for j in range(atoms):
            numeric_modes[cnt,0]= float(modes[i]['displacement'][j][0])
            numeric_modes[cnt,1]= float(modes[i]['displacement'][j][1])
            numeric_modes[cnt,2]= float(modes[i]['displacement'][j][2])
            cnt+=1
   
    return numeric_modes

if __name__ == "__main__":
    with open('../inputs.json') as f:
        inputs=json.load(f)
    if inputs["Molecule_data"]["Geom_flg"] == "PubChem":
        symbols,coordinates = get_geometry_pubchem(inputs["Molecule_data"]["Molecule"])
    elif inputs["Molecule_data"]["Geom_flg"] == "Initial":
        symbols,coordinates = get_geometry_file()
        generate_bondarr(symbols, coordinates)
    if inputs["Elec_structure"]["method"] == "QChem":
        import qchem as qc
        opt_coords, modes = qc.initial_conditions(symbols,coordinates,inputs['HPC']['cores'])
    elif inputs["Elec_structure"]["method"] == "PySCF":
        import elec_pyscf as pyscf
        pyscf.initial_conditions(symbols,coordinates)
    
   
    modes = modes['modes']
    num_modes = len(modes)
    natoms = int((num_modes+6)/3)
    modes=organise_modes(modes,natoms)
    masses = init.setup_masses(symbols)
    Px, Py, Pz = create_geom(natoms,num_modes,inputs["Molecule_data"]["Temp"],modes,masses,inputs["HPC"]["repeats"])
    momenta_checks(Px,Py,Pz,symbols,inputs["Molecule_data"]["Temp"])
    for j in range(inputs["HPC"]["repeats"]):
        with open('../repetitions/rep-'+str(j+1)+'/Geometry', 'w') as file:
            file.write(opt_coords)
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(natoms):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')
        
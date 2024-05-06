# File that handles outputs for all the different parts of the molecule: 
# 1. output xyz
# 2. output momentum
# 3. output forces and amplitudes
# 4. output json file for restarting

import init
import json


def output_molecule(molecule): 
    output_xyz(molecule)
    output_momenta(molecule)
    output_forces(molecule)
    molecule.to_json('../output/molecule.json')

def output_xyz(molecule):
    with open("../output/xyz.all", "a") as xyz_file:
        # Write the timestep
        xyz_file.write(f"Timestep: {molecule.timestep}\n")

        for atom_number, (symbol, (x, y, z)) in enumerate(zip(molecule.symbols, molecule.coordinates), start=1):
            xyz_file.write(f"{atom_number} {symbol}   {x:.15f}   {y:.15f}   {z:.15f}\n")

        # Write a dashed line as a separator
        xyz_file.write("-" * 40 + "\n")

def output_momenta(molecule):
    with open("../output/momenta.all", "a") as momenta_file:
        # Write the timestep
        momenta_file.write(f"Timestep: {molecule.timestep}\n")

        # Write the momenta
        momenta_file.write("Momenta:\n")
        for px, py, pz in molecule.momenta:
            momenta_file.write(f"{px:.15f}   {py:.15f}   {pz:.15f}\n")

        # Write a dashed line as a separator
        momenta_file.write("-" * 40 + "\n")

def output_forces(molecule):
    with open("../output/forces.all", "a") as forces_file:
        # Write the timestep
        forces_file.write(f"Timestep: {molecule.timestep}\n")

        # Write the amplitudes
        forces_file.write(f"Amplitudes: {molecule.amplitudes}\n")

        forces_file.write(f"Multiplicity: {molecule.multiplicity}\n")
        
        # Write the SCF energy
        forces_file.write(f"SCF Energy: {molecule.scf_energy}\n")

        # Write the forces
        forces_file.write("Forces:\n")
        for force in molecule.forces:
            forces_file.write(f"{force:.15f}\n")

        # Write dissociation information
        forces_file.write(f"Dissociation Flags: {molecule.dissociation_flags}\n")

        # Write a dashed line as a separator
        forces_file.write("-" * 40 + "\n")

def recombine_isolates(molecule_array):
    print(molecule_array[0].print_info())
    for i in range(1, len(molecule_array)):
        for j in range(0,len(molecule_array[i].symbols)):
 
            index = molecule_array[i].indexes[j] -1

            molecule_array[0].coordinates[index] =  molecule_array[i].coordinates[j]
            molecule_array[0].momenta[index] =  molecule_array[i].momenta[j]
            molecule_array[0].forces[index:index+3] =  molecule_array[i].forces[j:j+3]
    print('now: ')
    print(molecule_array[0].print_info())

    output_molecule(molecule_array[0])
            




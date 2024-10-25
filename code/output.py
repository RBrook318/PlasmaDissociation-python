"""
output.py

This module provides functionality for outputting the properties of a Molecule 
instance to various file formats. It includes functions to write the XYZ coordinates, 
momenta, forces, and additional attributes of the molecule to separate output files. 
The data is organized to facilitate analysis and visualization, with each output 
structured for easy parsing.

The following functions are available in this module:
- output_molecule: Aggregates output processes for a Molecule instance.
- output_xyz: Outputs the XYZ coordinates of the molecule to a file.
- output_momenta: Writes the momenta of the molecule to a file.
- output_forces: Outputs forces and additional properties of the molecule to a file.
"""

def output_molecule(molecule): 
    """
    Outputs the properties of a Molecule instance to various formats.

    This function calls the `output_xyz`, `output_momenta`, and `output_forces` 
    functions to write the coordinates, momenta, forces, and other properties of the 
    molecule to corresponding output files. It also serializes the Molecule instance 
    to a JSON file.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose properties are to be outputted.

    Functions
    ---------
    output_xyz - output.py
    output_momenta - output.py
    output_forces - output.py

    Returns
    -------
    None
    """
    output_xyz(molecule)
    output_momenta(molecule)
    output_forces(molecule)
    molecule.to_json('output/molecule.json')


def output_xyz(molecule):
    """
    Writes the XYZ coordinates of the Molecule instance to a file.

    This function appends the XYZ representation of the molecule's coordinates 
    along with the timestep to the "output/xyz.all" file. Each atom's information 
    is formatted for easy parsing.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose coordinates are to be written.

    Returns
    -------
    None
    """
    with open("output/xyz.all", "a") as xyz_file:
        # Write the timestep
        xyz_file.write(f"Timestep: {molecule.timestep}\n")

        for atom_number, (symbol, (x, y, z)) in enumerate(zip(molecule.symbols, molecule.coordinates), start=1):
            xyz_file.write(f"{atom_number} {symbol}   {x:.15f}   {y:.15f}   {z:.15f}\n")

        # Write a dashed line as a separator
        xyz_file.write("-" * 40 + "\n")


def output_momenta(molecule):
    """
    Writes the momenta of the Molecule instance to a file.

    This function appends the momenta of the molecule to the "output/momenta.all" 
    file, along with the current timestep.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose momenta are to be written.

    Returns
    -------
    None
    """
    with open("output/momenta.all", "a") as momenta_file:
        # Write the timestep
        momenta_file.write(f"Timestep: {molecule.timestep}\n")

        # Write the momenta
        momenta_file.write("Momenta:\n")
        for px, py, pz in molecule.momenta:
            momenta_file.write(f"{px:.15f}   {py:.15f}   {pz:.15f}\n")

        # Write a dashed line as a separator
        momenta_file.write("-" * 40 + "\n")


def output_forces(molecule):
    """
    Writes the forces and additional properties of the Molecule instance to a file.

    This function appends the forces, amplitudes, multiplicity, SCF energy, 
    and dissociation flags of the molecule to the "output/forces.all" file, 
    along with the current timestep.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose forces and properties are to be written.

    Returns
    -------
    None
    """
    with open("output/forces.all", "a") as forces_file:
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

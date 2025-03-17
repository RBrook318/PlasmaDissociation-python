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
import os
import numpy as np
import matplotlib.pyplot as plt

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
    output_time(molecule)
    molecule.to_json('output/molecule.json')

def run_checks(molecule):
    """
    Runs the overall checks for the molecule.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose coordinates are to be written.
    
    
    """
    forces_magnitudes(molecule)

def forces_magnitudes(molecule):
    # List to store indices of atoms that have not dissociated
    shrunk_index = []
    natoms = len(molecule.symbols)
    
    # Identifying atoms that have not dissociated
    for i in range(natoms): 
        if molecule.dissociation_flags[i] == 'NO':
            shrunk_index.append(i)
    
    # Loop through atoms and calculate the force magnitudes
    for j in range(len(shrunk_index)):
        # Calculate the force magnitude (squared sum of force components)
        force_magnitude = np.sqrt(np.sum(molecule.forces[j, :3,:]**2))  
        print(type(force_magnitude))
        # Save the force magnitude to a file named after the atom index
        with open(f"checks/forcemagnitude_{shrunk_index[j]+1}.out", "a") as out_file:
            out_file.write(f"{force_magnitude}\n")
    
    # After all files are written, we can plot the force magnitude evolution
    plot_force_magnitudes()

def plot_force_magnitudes(directory="checks"):
    # Attempt to read fragment-time.out if it exists
    fragment_timesteps = []
    fragment_file = os.path.join(directory, "fragment-time.out")
    
    if os.path.exists(fragment_file):
        with open(fragment_file, 'r') as f:
            for line in f:
                # Convert each timestep to an integer and store it
                fragment_timesteps.extend(map(lambda x: int(float(x)), line.split()))

    
    # Loop over all relevant files in the directory
    for filename in os.listdir(directory):
        if filename.startswith("forcemagnitude_") and filename.endswith(".out"):
            # Extract atom index from filename
            atom_index = int(filename.split("_")[1].split(".")[0])

            # Read the force magnitudes from the file
            force_magnitudes = []
            with open(os.path.join(directory, filename), 'r') as f:
                for line in f:
                    # Convert each line (force magnitude) to a float and add it to the list
                    force_magnitudes.append(float(line.strip()))
            # Create an x-axis representing the position in the file (time steps)
            time_steps = np.arange(len(force_magnitudes))

            # Plot the force magnitudes over time
            plt.figure(figsize=(8, 6))
            plt.plot(time_steps, force_magnitudes, label=f'Atom {atom_index}', color='b')
            plt.xlabel('Time Step (Position in File)')
            plt.ylabel('Force Magnitude')
            plt.title(f'Force Magnitude Evolution for Atom {atom_index}')
            plt.legend()
            plt.grid(True)

            # If fragment timesteps exist, mark them on the graph
            if fragment_timesteps:
                for t in fragment_timesteps:
                    if 0 <= t < len(force_magnitudes):  # Ensure the timestep is within range
                        plt.axvline(x=t, color='r', linestyle='dashed', alpha=0.7, label="Fragmentation Event")

                # Ensure the legend only contains one label for fragmentation events
                handles, labels = plt.gca().get_legend_handles_labels()
                unique_labels = dict(zip(labels, handles))
                plt.legend(unique_labels.values(), unique_labels.keys())

            # Save the plot as an image
            plt.savefig(f"{directory}/force_magnitude_atom_{atom_index}.png")
            plt.close()
    
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
        # List to store indices of atoms that have not dissociated
    shrunk_index = []
    natoms = len(molecule.symbols)
    
    # Identifying atoms that have not dissociated
    for i in range(natoms): 
        if molecule.dissociation_flags[i] == 'NO':
            shrunk_index.append(i)
    with open("output/forces.all", "a") as forces_file:
        # Write the timestep
        forces_file.write(f"Timestep: {molecule.timestep}\n")

        # Write the amplitudes
        forces_file.write(f"Amplitudes: {molecule.amplitudes}\n")

        forces_file.write(f"Multiplicity: {molecule.multiplicity}\n")
        
        # Write the SCF energy
        forces_file.write(f"SCF Energy: {molecule.scf_energy}\n")

        forces_file.write("Forces:\n")
        for atom_idx in range(molecule.forces.shape[0]):  # Loop over atoms
            forces_file.write(f"Atom {shrunk_index[atom_idx]+ 1}:\n")
            for state_idx in range(molecule.forces.shape[2]):  # Loop over states
                fx, fy, fz = molecule.forces[atom_idx, :, state_idx]
                forces_file.write(f"State {state_idx + 1}: {fx:.8f} {fy:.8f} {fz:.8f}\n")
            # Check if the coupling array is all zeros
            if not np.all(molecule.coupling == 0):  # Only write coupling if nonzero
                cx, cy, cz = molecule.coupling[atom_idx, :]
                forces_file.write(f"Coupling: {cx:.8f} {cy:.8f} {cz:.8f}\n")

        # Write dissociation information
        forces_file.write(f"Dissociation Flags: {molecule.dissociation_flags}\n")

        # Write a dashed line as a separator
        forces_file.write("-" * 40 + "\n")

def output_time(molecule):
     with open("output/time.all", "a") as times_file:
        times_file.write(f"Timestep: {molecule.timestep}\n")
        times_file.write(f"Time of this step: {molecule.time[1]}\n")
        times_file.write(f"Time of QChem this step: {molecule.time[0]}\n")
        times_file.write(f"Time of python this step: {molecule.time[1]-molecule.time[0]}\n")
        times_file.write(f"Total time of propagation: {molecule.time[3]}s (average: {molecule.time[3]/molecule.time[4]}s) \n")
        times_file.write(f"Total time of QChem: {molecule.time[2]} (average: {molecule.time[2]/molecule.time[4]}s)\n")
        times_file.write("-" * 40 + "\n")

def output_fragment_time(timestep):
    with open("checks/fragment-time.out", "a") as time_file:
        time_file.write(f"{timestep}\n")
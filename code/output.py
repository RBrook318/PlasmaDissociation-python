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
import time
import global_vars as gv

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

def run_checks(molecule, old_coordinates):
    """
    Runs the overall checks for the molecule.

    Parameters
    ----------
    molecule : Molecule
        The Molecule instance whose coordinates are to be written.
    old_coordinates : array-like
        The previous coordinates of the molecule.
    """

    
    # Timing forces_magnitudes

    forces_magnitudes(molecule)

   
    
    # Timing check_total_energy

    check_total_energy(molecule)

    
    # Timing displacement_vector
    displacement_vector(molecule, old_coordinates)
    
    end_time = time.time()

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
    natoms = len(molecule.symbols)
    
    # Identifying atoms that have not dissociated
    shrunk_index = np.where(np.array(molecule.dissociation_flags) == 'NO')[0]

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
            if len(molecule.scf_energy) > 1: # Only write coupling if nonzero
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
        times_file.write(f"Total time of propagation: {molecule.time[3]}s (average: {molecule.time[3]/molecule.timestep/gv.timestep}s) \n")
        times_file.write(f"Total time of QChem: {molecule.time[2]} (average: {molecule.time[2]/(molecule.timestep/gv.timestep)}s)\n")
        times_file.write("-" * 40 + "\n")

def output_fragment_time(timestep):
    with open("checks/fragment-time.out", "a") as time_file:
        time_file.write(f"{timestep}\n")

def check_total_energy(molecule):
    # Fixed directory path
    energy_dir = "checks/energy"
    
    # Ensure the directory exists
    os.makedirs(energy_dir, exist_ok=True)

    # Identifying non-dissociated atoms
    shrunk_index = np.where(np.array(molecule.dissociation_flags) == 'NO')[0]
    # Compute potential energy
    potential_energy = np.sum(molecule.scf_energy * np.abs(molecule.amplitudes) ** 2)

    # Compute kinetic energy (vectorized)
    kinetic_energy = 0.5 * np.sum((molecule.momenta[shrunk_index] ** 2) / molecule.masses[shrunk_index, np.newaxis])

    # Total energy
    total_energy = potential_energy + kinetic_energy

    # Write total energy to file (append mode to track evolution)
    energy_file = os.path.join(energy_dir, "total_energy.out")
    with open(energy_file, "a") as f:
        f.write(f"{total_energy:.10f}\n")

    # Read total energy values
    total_energy_values = []
    with open(energy_file, "r") as f:
        for line in f:
            total_energy_values.append(float(line.strip()))

    # Read fragmentation timesteps
    fragment_timesteps = []
    fragment_file = "checks/fragment-time.out"  # Fixed location
    
    if os.path.exists(fragment_file):
        with open(fragment_file, 'r') as f:
            for line in f:
                fragment_timesteps.extend(map(lambda x: int(float(x)), line.split()))

    # Generate time steps based on the energy values read
    time_steps = np.arange(len(total_energy_values))

    # Plot change in total energy from the first entry
    plt.figure(figsize=(8, 6))
    
    if len(total_energy_values) > 1:
        initial_energy = total_energy_values[0]
        delta_energy_values = [value - initial_energy for value in total_energy_values]
        plt.plot(time_steps, delta_energy_values, label="Change in Total Energy", color='b')
        plt.xlabel('Time Step')
        plt.ylabel('Change in Total Energy')
        plt.title('Change in Total Energy Evolution')
        plt.legend()
        plt.grid(True)

        # Mark fragmentation events
        if fragment_timesteps:
            for t in fragment_timesteps:
                if 0 <= t < len(delta_energy_values):  # Ensure the timestep is within range
                    plt.axvline(x=t, color='r', linestyle='dashed', alpha=0.7, label="Fragmentation Event")

            # Ensure the legend only contains one label for fragmentation events
            handles, labels = plt.gca().get_legend_handles_labels()
            unique_labels = dict(zip(labels, handles))
            plt.legend(unique_labels.values(), unique_labels.keys())
        
        # Save the plot
        plt.savefig(os.path.join(energy_dir, "change_in_total_energy_evolution.png"))
        plt.close()
    
    else:
        print("Not enough data to plot change in total energy.")
    
def forces_magnitudes(molecule):
    # Define the fixed directory path
    forces_dir = "checks/forces"
    
    # Ensure the directory exists
    os.makedirs(forces_dir, exist_ok=True)

    # List to store indices of atoms that have not dissociated
 
    shrunk_index = np.where(np.array(molecule.dissociation_flags) == 'NO')[0]

    # Loop through non-dissociated atoms and calculate force magnitudes
    for j in range(len(shrunk_index)):
        # Calculate force magnitude (L2 norm)
        force_magnitude = np.sqrt(np.sum(molecule.forces[j, :3, :]**2))

        # Save the force magnitude to a file named after the atom index
        with open(os.path.join(forces_dir, f"forcemagnitude_{shrunk_index[j]+1}.out"), "a") as out_file:
            out_file.write(f"{force_magnitude}\n")

    # Plot force magnitudes after writing files
    if gv.checks > 1:
        plot_force_magnitudes(forces_dir)

def plot_force_magnitudes(directory="checks/forces"):

    # Ensure the directory exists
    os.makedirs(directory, exist_ok=True)

    # Read fragment-time.out if it exists
    fragment_timesteps = []
    fragment_file = os.path.join(directory, "../fragment-time.out")  # Assumed to be in 'checks/'
    
    if os.path.exists(fragment_file):
        with open(fragment_file, 'r') as f:
            for line in f:
                fragment_timesteps.extend(map(lambda x: int(float(x)), line.split()))

    # Loop over all relevant files in the directory
    for filename in os.listdir(directory):
        if filename.startswith("forcemagnitude_") and filename.endswith(".out"):
            # Extract atom index from filename
            atom_index = int(filename.split("_")[1].split(".")[0])

            # Read the force magnitudes from the file
            force_magnitudes = np.loadtxt(os.path.join(directory, filename))
            force_magnitudes = np.atleast_1d(force_magnitudes)
            # Create an x-axis representing time steps
            time_steps = np.arange(len(force_magnitudes))

            # Plot force magnitudes over time
            plt.figure(figsize=(8, 6))
            plt.plot(time_steps, force_magnitudes, label=f'Atom {atom_index}', color='b')
            plt.xlabel('Time Step')
            plt.ylabel('Force Magnitude')
            plt.title(f'Force Magnitude Evolution for Atom {atom_index}')
            plt.legend()
            plt.grid(True)

            # Mark fragmentation events
            if fragment_timesteps:
                for t in fragment_timesteps:
                    if 0 <= t < len(force_magnitudes):  
                        plt.axvline(x=t, color='r', linestyle='dashed', alpha=0.7, label="Fragmentation Event")

                # Ensure the legend contains unique labels
                handles, labels = plt.gca().get_legend_handles_labels()
                unique_labels = dict(zip(labels, handles))
                plt.legend(unique_labels.values(), unique_labels.keys())

            # Save the plot
            plt.savefig(os.path.join(directory, f"force_magnitude_atom_{atom_index}.png"))
            plt.close()

def displacement_vector(molecule, old_coordinates):
    # Define the directory
    disp_dir = "checks/displacements"
    os.makedirs(disp_dir, exist_ok=True)

    # Identify non-dissociated atoms
    shrunk_index = np.where(np.array(molecule.dissociation_flags) == 'NO')[0]

    # Compute displacement for each atom
    displacements = molecule.coordinates - old_coordinates

    displacement_magnitudes = np.sqrt(np.sum(displacements ** 2, axis=1))  # Per-atom magnitude

    # Compute total displacement
    total_displacement = np.sum(displacement_magnitudes)

    # Save individual displacements
    for i in shrunk_index:
        with open(os.path.join(disp_dir, f"displacement_{i+1}.out"), "a") as out_file:
            out_file.write(f"{displacement_magnitudes[i]:.10f}\n")    



    # Save total displacement
    with open(os.path.join(disp_dir, "total_displacement.out"), "a") as out_file:
        out_file.write(f"{total_displacement:.10f}\n")


     

    # Plot displacement evolution
    if gv.checks > 1:
        plot_displacements(disp_dir)

def plot_displacements(directory="checks/displacements"):
    # Ensure directory exists
    os.makedirs(directory, exist_ok=True)

    # Load fragmentation events if they exist
    fragment_timesteps = []
    fragment_file = os.path.join(directory, "../fragment-time.out")  # Assuming stored in "checks/"
    if os.path.exists(fragment_file):
        with open(fragment_file, 'r') as f:
            fragment_timesteps = [int(float(x)) for x in f.read().split()]

    # Plot individual atom displacements
    for filename in os.listdir(directory):
        if filename.startswith("displacement_") and filename.endswith(".out"):
            atom_index = int(filename.split("_")[1].split(".")[0])
            displacement_data = np.loadtxt(os.path.join(directory, filename))
            displacement_data = np.atleast_1d(displacement_data) 
            time_steps = np.arange(len(displacement_data))

            plt.figure(figsize=(8, 6))
            plt.plot(time_steps, displacement_data, label=f'Atom {atom_index}', color='b')
            plt.xlabel('Time Step')
            plt.ylabel('Displacement Magnitude')
            plt.title(f'Displacement Evolution for Atom {atom_index}')
            plt.legend()
            plt.grid(True)

            # Mark fragmentation events
            if fragment_timesteps:
                for t in fragment_timesteps:
                    if 0 <= t < len(displacement_data):
                        plt.axvline(x=t, color='r', linestyle='dashed', alpha=0.7, label="Fragmentation Event")

                # Ensure only one fragmentation event label appears
                handles, labels = plt.gca().get_legend_handles_labels()
                unique_labels = dict(zip(labels, handles))
                plt.legend(unique_labels.values(), unique_labels.keys())

            plt.savefig(os.path.join(directory, f"displacement_atom_{atom_index}.png"))
            plt.close()

    # Plot total displacement
    total_disp_file = os.path.join(directory, "total_displacement.out")
    if os.path.exists(total_disp_file):
        total_displacements = np.loadtxt(total_disp_file)
        total_displacements = np.atleast_1d(total_displacements) 
        time_steps = np.arange(len(total_displacements))

        plt.figure(figsize=(8, 6))
        plt.plot(time_steps, total_displacements, label="Total Displacement", color='g')
        plt.xlabel("Time Step")
        plt.ylabel("Total Displacement")
        plt.title("Total Displacement Over Time")
        plt.legend()
        plt.grid(True)

        # Mark fragmentation events
        if fragment_timesteps:
            for t in fragment_timesteps:
                if 0 <= t < len(total_displacements):
                    plt.axvline(x=t, color='r', linestyle='dashed', alpha=0.7, label="Fragmentation Event")

            handles, labels = plt.gca().get_legend_handles_labels()
            unique_labels = dict(zip(labels, handles))
            plt.legend(unique_labels.values(), unique_labels.keys())

        plt.savefig(os.path.join(directory, "total_displacement.png"))
        plt.close()
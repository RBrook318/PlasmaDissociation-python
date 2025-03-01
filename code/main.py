"""
main.py

This script serves as the main entry point for running molecular simulations using the 
custom-built framework. It orchestrates the entire simulation process, from initializing 
the molecular structure to processing results. 

Key Functionalities:
--------------------
- Loads input parameters from a JSON configuration file.
- Initializes the molecular structure based on specified simulation parameters.
- Runs electronic structure calculations and properties evaluations.
- Handles restart functionality for previously interrupted simulations.
- Outputs molecular data at each timestep for analysis.
- Measures and records the execution time of the simulation.
- Processes and collates the results of the simulation at the end.

Input:
------
- A JSON file (`../inputs.json`) that contains configuration settings for the simulation.
  This file must specify the following parameters:
  - `repeats`: Number of repetitions for the simulation.
  - `cores`: Number of CPU cores to use for calculations.
  - `States`: Number of states to consider in the simulation.
  - `Branches`: Whether to allow for cloning ().
  - `Multiplicity`: Initial spin multiplicity for the molecule.
  - `Timestep`: Time increment for the simulation steps.
  - `Tot_timesteps`: Total number of timesteps to run.
  - `Geom_start`: If reading from a folder of geometries, which to start from.
  - `Spin_flip`: Boolean to indicate if spin-flip is to be used.
  - `method`: The electronic structure package to use.

Output:
-------
- Outputs molecular data to the `output/` directory, including:
  - `xyz.all`: A file containing the coordinates trajectory data.
  - `moementa.all`: A file containing the momenta trajectory data.
  - `forces.all`: A file containing the forces, multiplicity, dissociation flag of the trajectory.
  - `time.out`: A file recording the execution time of the simulation.
  - JSON files representing the molecular structure at each timestep.

Usage:
------
Run the script from the command line. Ensure that the input JSON file is correctly formatted 
and accessible at the specified path. The script will handle initialization, execution, 
and result processing, with proper checks for whether a simulation can be restarted.

Exceptions:
-----------
- The script checks if a simulation has already completed by examining the `xyz.all` file.
  If the run is complete, it will exit and not process results to avoid unnecessary computation.
- It raises a `ValueError` if any expected files or directories are missing.

Author: Ryan Brook
Date: 30/10/2024
"""

import os
import json
import time
from init import Molecule, initialize_structure, create_empty_molecule
import elec
import prop
import output as out
import result

# Constants for file paths
INPUTS_FILE = '../inputs.json'
XYZ_FILE = 'output/xyz.all'
MOLECULE_JSON_FILE = 'output/molecule.json'
TIME_OUTPUT_FILE = 'output/time.out'


def load_inputs():
    """Load simulation parameters from a JSON file."""
    with open(INPUTS_FILE) as f:
        return json.load(f)

def check_restart(endstep, increment):
    """Check if the simulation can be restarted based on the existing XYZ file."""
    if os.path.exists(XYZ_FILE):
        with open(XYZ_FILE, "r") as xyz_file:
            lines = xyz_file.readlines()
        for line in reversed(lines):
            if "Timestep:" in line:
                break
        line = line.split()
        third_last_line = int(line[1]) / increment
        print("Time step: ", third_last_line)
        return 'YES' if third_last_line < endstep else 'NO'
    return 'NO'

def initialize_simulation(inputs, restart,basis):
    """Initialize the molecular structure and handle restart logic."""
    nstates = inputs["run"]["States"]
    spin_flip = inputs["run"]["Spin_flip"]
    mult = inputs["run"]["Multiplicity"]    
    time1 =time.time()
    if restart == 'NO':
        molecule1 = initialize_structure(nstates, spin_flip, mult)
        n = len(molecule1.symbols)
        molecule2 = create_empty_molecule(n, nstates, spin_flip)
        molecule1 = elec.run_elec_structure(molecule1, inputs["setup"]["cores"], nstates, spin_flip, inputs["run"]["method"], Guess=False, basis = basis)
        time2= time.time()
        molecule1.time[1] = time2-time1
        molecule1.time[3] += time2-time1
        molecule1.time[4] = molecule1.timestep/inputs["run"]["Timestep"]
        out.output_molecule(molecule1)
        molecule1.time[0] = 0
        return molecule1, molecule2, 1, True
    # Restart logic
    if os.path.exists(MOLECULE_JSON_FILE):
        molecule1 = Molecule.from_json(MOLECULE_JSON_FILE)
        n = len(molecule1.symbols)
        molecule2 = create_empty_molecule(n, nstates, spin_flip)
        startstep = molecule1.timestep / inputs["run"]["Timestep"]
        return molecule1, molecule2, startstep, False
    
    # If JSON file doesn't exist, initialize a new structure
    molecule1 = initialize_structure(nstates, spin_flip, mult)
    n = len(molecule1.symbols)
    molecule2 = create_empty_molecule(n, nstates, spin_flip)
    molecule1 = elec.run_elec_structure(molecule1, inputs["setup"]["cores"], nstates, spin_flip, inputs["run"]["method"], Guess=False, basis = basis)
    time2= time.time()
    molecule1.time[1] = time2-time1
    molecule1.time[3] += time2-time1
    molecule1.time[4] = molecule1.timestep/inputs["run"]["Timestep"]
    out.output_molecule(molecule1)
    molecule1.time[0] = 0
    return molecule1, molecule2, 1, True

def run_simulation(inputs, startstep, endstep, molecule1, molecule2, guess, basis):
    """Run the main simulation loop for the specified number of timesteps."""
    n = len(molecule1.symbols)
    for i in range(int(startstep), endstep + 1):
        time1 = time.time()
        molecule2 = prop.prop_1(molecule1, molecule2, n, inputs["run"]["States"], inputs["run"]["Timestep"])
        molecule2 = elec.run_elec_structure(molecule2, inputs["setup"]["cores"], inputs["run"]["States"], inputs["run"]["Spin_flip"], inputs["run"]["method"], Guess=guess,basis=basis)
        molecule2.coupling = prop.align_coupling(molecule1.coupling, molecule2.coupling)
        molecule1.time[0] = molecule2.time[0]
        molecule1.time[2] = molecule2.time[2]
        molecule1.elecinfo = molecule2.elecinfo

        molecule1 = prop.prop_2(molecule1, molecule2, n, inputs["run"]["States"], inputs["run"]["Timestep"])
        molecule1, dissociated = prop.fragments(molecule1, inputs["run"]["Spin_flip"])
        molecule1 = prop.prop_diss(molecule1, inputs["run"]["Timestep"])
        time2= time.time()
        molecule1.time[1] = time2-time1
        molecule1.time[3] += time2-time1
        molecule1.time[4] = molecule1.timestep/inputs["run"]["Timestep"]
        out.output_molecule(molecule1)
        out.run_checks(molecule1)
        molecule1.time[0] = 0
        guess = dissociated == 0  # Update guess based on dissociation
        
def main():

    inputs = load_inputs()

    # Check basic arguments
    basis = inputs["run"]["Basis"]
    endstep = inputs["run"]["Tot_timesteps"]
    increment = inputs["run"]["Timestep"]
    restart = check_restart(endstep, increment)

    molecule1, molecule2, startstep, guess = initialize_simulation(inputs, restart,basis)

    run_simulation(inputs, startstep, endstep, molecule1, molecule2, guess, basis)

    result.process_results()

if __name__ == "__main__":
    main()



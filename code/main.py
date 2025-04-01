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
from init import Molecule, initialize_structure, create_empty_molecule, initialize_simulation
import elec
import prop
import output as out
import result
import global_vars as gv

# def initialize_simulation(restart):
#     """Initialize the molecular structure and handle restart logic."""
   
#     time1 =time.time()
#     if restart == 'NO':
#         time2= time.time()
#         molecule1.time[1] = time2-time1
#         molecule1.time[3] += time2-time1
#         molecule1.time[4] = molecule1.timestep/gv.timestep
        
#         molecule1.time[0] = 0
#         return molecule1, molecule2, 1, True
    
#     # If JSON file doesn't exist, initialize a new structure
#     molecule1 = initialize_structure()
#     n = len(molecule1.symbols)
#     molecule2 = create_empty_molecule(n)
#     molecule1 = elec.run_elec_structure(molecule1, Guess=False)
#     time2= time.time()
#     molecule1.time[1] = time2-time1
#     molecule1.time[3] += time2-time1
#     molecule1.time[4] = molecule1.timestep/gv.timestep
#     out.output_molecule(molecule1)
#     molecule1.time[0] = 0
#     return molecule1, molecule2, 1, True

def run_simulation(startstep, molecule1, molecule2):
    """Run the main simulation loop for the specified number of timesteps."""
    n = len(molecule1.symbols)
    guess = True

    for i in range(int(startstep)-1, gv.tot_timesteps + 1):
        time1 = time.time()
        old_coordinates = molecule1.coordinates.copy()
        print(molecule1.coordinates)
        molecule2 = prop.prop_1(molecule1, molecule2)
        
        molecule2 = elec.run_elec_structure(molecule2, Guess=guess)

        molecule1.time[0] = molecule2.time[0]
        molecule1.time[2] = molecule2.time[2]

        molecule1.elecinfo = molecule2.elecinfo
        prop1t = time.time()
        molecule1 = prop.prop_2(molecule1, molecule2)
        prop2t = time.time()
        print("Prop2 time = ", prop2t-prop1t)
        remove1t = time.time()
        if gv.remove_atoms == 1:
            molecule1, dissociated = prop.fragments(molecule1)
            molecule1 = prop.prop_diss(molecule1)
            guess = dissociated == 0  # Update guess based on dissociation
        remove2t = time.time()
        print("Remove time = ", remove2t - remove1t)
        
        check1t = time.time()
        if gv.checks > 0:
            out.run_checks(molecule1,old_coordinates)
        time2= time.time()
        molecule1.time[1] = time2-time1
        molecule1.time[3] += time2-time1
        check2t = time.time()
        print("Check time = ", check2t-check1t)
        out.output_molecule(molecule1)

        molecule1.time[0] = 0

def main():


    gv.load_global_variables()
    # Check basic arguments

    # restart = init.check_restart(gv.tot_timesteps, gv.timestep)

    molecule1, molecule2, startstep = initialize_simulation()

    out.output_molecule(molecule1)
    
    run_simulation(startstep,molecule1, molecule2)

    result.process_results()

if __name__ == "__main__":
    main()



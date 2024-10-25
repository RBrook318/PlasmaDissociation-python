# Plasma Dissociation Code

## Description

This code runs a molecular dynamics simulation for molecules on (currently only) 1 potential energy surface. 

Initial geometry will be optimised, and initial momentum generated from a Boltzmann distribution and then propagated via a Verlett Algorithm. 

The acceleration for these propagations can be found via the electronic forces (F=ma) acting on each atom within the molecule as calculated by electronic structure packages (QChem/PySCF).

At the end of the dynamics, the trajectory is sacnned through to find the time at which any bonds were broken and automatically generates plot to display the average number of bonds broken per trajectory. 


## Installation
If not already made create a new conda environment for the project using the command 

conda create --name <project environement name>

Then following packages can then be added to the environment to ensure proper compilation.
Python package requirements
- numpy -- pip install pyinstaller
- scipy -- pip install scipy
- pubchempy -- pip install pubchempy
- subrprocess
- re 
- math
- sys
- socket
- json
- shutil
- networkx
- pyqchem -- pip install pyqchem (if using QChem)
- PySCF -- pip install PySCF (if using PySCF)

# For initial installation
pip install numpy scipy pubchempy subprocess re math sys socket json shutil networkx
# For QChem users
pip install pyqchem

Note that when I installed pyqchem the tools folder of packages did not install correctly. I have added all the relevant files in the folder tools. This needs to put into the folder .conda/env/$project_environment_name$/lib/python3.#/site-packages/pyqchem

# For PySCF users
pip install pyscf




# Running on university of Leeds HPC only 
I have also found it useful to add a shortcut in my .bashrc file to quickly load the appropriate modules on arc. The foloowing can be added 

scatter(){
	module add qchem
	module swith gnu intel
	module add intel 
	module switch intel/18.0.2
	module add mkl
	module add anaconda
	source activate scatter
}

The conda environment I have named scatter and this is also the shortcut command to load the required modules. 

## Usage
To run the program the mkl and qchem module needs to be loaded. The run.py script will set up all the relevant folders submitting an initial setup job and the subsequent run.


Here is an example of the `input.json` file required to run the program:

- **`setup` section**: Configures the main directory, core usage, and run options.
- **`run` section**: Defines molecular dynamics parameters, including temperature, number of states, and method.

```json
{
    "setup": {
        "Runfolder": "methane-test",     // Name of the main run folder
        "repeats": 1,                    // Number of repeated simulations
        "cores": 8,                      // Number of CPU cores to use
        "restart": "NO",                 // Set to "YES" if restarting a run
        "initre": 0                      // Initialization step (0 = default start)
    },
    "run": {
        "Geom_flg": 0,                   // Geometry flag (0 = take from folder, 1 = create new (setup job))
        "Molecule": "Methane",           // Name of the molecule (name of folder, or molecule ID for making new geometries)
        "States": 1,                     // Number of electronic states to simulate
        "Multiplicity": 3,               // Initial spin multiplicity of the molecule
        "Temp": 1000,                    // Temperature in Kelvin
        "Cutoff": 1,                     // Cutoff for highest percentage of energies taken as 
        "Branches": 0,                   // Number of branches for dynamics (0 = none)
        "Timestep": 10,                  // Time step for MD simulation (in a.t.u)
        "Tot_timesteps": 2100,           // Total number of time steps (2100 a.t.u = 500 fs)
        "Geom_start": 1,                 // Starting geometry (1 = optimized)
        "Spin_flip": 0,                  // Whether spin flip DFT is enabled (0 = off)
        "method": "PySCF",               // Method for electronic structure calculation (PySCF or QChem)
        "GPU": 0                         // Use GPU if available (0 = no, 1 = yes)
    }
}
## Contributing

[Add guidelines for contributing to your project]

## License

[Add information about the license for your project]

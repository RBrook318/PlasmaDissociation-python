# Plasma Dissociation Code

## Description

[Add a brief description of your project here]

## Installation
If not already made create a new conda environment for the project using the command 

conda create --name <project environement name>

Then following packages can then be added to the environment to ensure proper compilation.
Python package requirements
- numpy -- pip install pyinstaller
- scipy -- pip install scipy
- subrprocess
- re 
- math
- sys
- socket
- json
- shutil
- pyinstaller -- pip install pyinstaller
- chemformula -- pip install chemformula 
- pyqchem -- pip install pyqchem 

Note that when I installed pyqchem the tools folder of packages did not install correctly. I have added all the relevant files in the folder tools. This needs to put into the folder .conda/env/$project_environment_name$/lib/python3.#/site-packages/pyqchem

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

The conda environment I have named scatter and this is also the shortcut command to load the required modules

## Usage
To run the program the mkl and qchem module needs to be loaded. The run.py script will set up all the relevant folders submitting an initial setup job and the subsequent run.

## Contributing

[Add guidelines for contributing to your project]

## License

[Add information about the license for your project]

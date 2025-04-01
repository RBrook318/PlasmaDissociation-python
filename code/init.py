"""
init.py

This module provides utilities for initializing molecular dynamics simulation structures, including 
classes and functions for handling molecular properties such as symbols, coordinates, forces, momenta, 
and more. 

Classes:
    NumpyEncoder: Custom JSON encoder for handling numpy float and complex types.
    Molecule: Represents a molecule with attributes such as symbols, coordinates, forces, and more.

Functions:
    create_empty_molecule: Creates a molecule instance with default attributes.
    initialize_structure: Reads molecular structure data from a file and initializes a molecule instance.
    create_molecule: Wrapper function for molecule creation with optional structure reading.
    setup_masses: Returns atomic masses for a given set of symbols in atomic mass units (amu).
"""

import numpy as np
import json
import os
np.set_printoptions(precision=30)
import global_vars as gv
import elec

class NumpyEncoder(json.JSONEncoder):
    """
    A JSON encoder to handle numpy float and complex types for JSON serialization.

    Methods
    -------
    default(obj):
        Overrides default method to check for numpy float and complex types and convert them to standard
        Python float or complex representations.
    """
    def default(self, obj):
        if isinstance(obj, np.float64):
            return float(obj)
        return super(NumpyEncoder, self).default(obj)

class Molecule:
    """
    Represents a molecule with its properties and methods to manipulate them.

    Attributes
    ----------
    symbols : list[str]
        Chemical symbols for each atom.
    coordinates : np.ndarray
        Array of atomic positions in Cartesian coordinates.
    momenta : np.ndarray, optional
        Array of atomic momenta, if provided.
    scf_energy : np.ndarray
        Array of SCF energy levels.
    forces : np.ndarray, optional
        Array of forces acting on atoms, if provided.
    amplitudes : np.ndarray
        Array of complex amplitudes representing quantum state amplitudes.
    timestep : int
        The current time step in the simulation.
    multiplicity : int
        The molecule's spin multiplicity.
    dissociation_flags : list[str]
        Flags indicating the dissociation status of each atom.
    elecinfo : any, optional
        Additional electronic information.
    masses : np.ndarray
        Atomic masses of the elements.

    Methods
    -------
    update_*(new_info):
        Updates * of the molecule class with new_info.
    print_info():
        Prints all molecular attributes for inspection.
    copy():
        Returns a deep copy of the molecule instance.
    to_dict():
        Converts the molecule attributes to a dictionary format.
    to_json(filename):
        Saves the molecule attributes as a JSON file.
    from_dict(data):
        Class method to initialize a Molecule instance from a dictionary.
    from_json(filename):
        Class method to initialize a Molecule instance from a JSON file.
    """
    def __init__(self, symbols, coordinates, momenta=None, scf_energy=None, forces=None, amplitudes=None, timestep=0, multiplicity=5, dissociation_flags=None, elecinfo=None, masses=None, coupling=None,time=None):
        self.symbols = symbols
        self.coordinates = np.array(coordinates, dtype=np.float64)
        self.momenta = np.array(momenta, dtype=np.float64) if momenta is not None else None
        self.scf_energy = np.array(scf_energy, dtype=np.float64)
        self.forces = np.array(forces, dtype=np.float64) if forces is not None else None
        self.amplitudes = np.array(amplitudes, dtype=np.complex256) if amplitudes is not None else np.array([1.0 + 0.0j], dtype=np.complex256)
        self.timestep = timestep
        self.multiplicity = multiplicity
        self.dissociation_flags = dissociation_flags or ["NO"] * len(symbols)
        self.elecinfo = elecinfo if elecinfo is not None else None
        self.masses = np.array(masses, dtype=np.float64) if masses is not None else np.zeros(len(symbols), dtype=np.float64)
        self.coupling = np.array(coupling, dtype=np.float64) if coupling is not None else None
        self.time = np.array(time) if time is not None else None

    def update_symbols(self, new_symbols):
        self.symbols = new_symbols

    def update_coordinates(self, new_coordinates):
        self.coordinates = np.array(new_coordinates)

    def update_momenta(self, new_momenta):
        self.momenta = np.array(new_momenta)

    def update_scf_energy(self, new_scf_energy):
        self.scf_energy = np.array(new_scf_energy)

    def update_forces(self, new_forces):
        self.forces = np.array(new_forces)

    def update_amplitudes(self, new_amplitudes):
        self.amplitudes = np.array(new_amplitudes)

    def update_timestep(self, new_timestep):
        self.timestep = new_timestep

    def update_multiplicity(self, multiplicity):
        self.multiplicity = multiplicity
    
    def update_dissociation_flags(self, new_flags):
        self.dissociation_flags = new_flags

    def update_elecinfo(self, new_elecinfo):
        self.elecinfo = new_elecinfo

    def update_masses(self, new_masses):
        self.masses = np.array(new_masses, dtype=np.float64)
    
    def update_coupling(self, new_coupling):
        self.coupling = np.array(new_coupling, dtype=np.float64)
    

    def print_info(self):
        print("Symbols:", self.symbols)
        print("Coordinates:")
        print(self.coordinates)
        print("Momenta:")
        print(self.momenta)
        print("SCF Energy:", self.scf_energy)
        print("Forces:")
        print(self.forces)
        print("Amplitudes:")
        print(self.amplitudes)
        print("Timestep:", self.timestep)
        print("Multiplicity:", self.multiplicity)
        print("Dissociation Flags:", self.dissociation_flags)
        print("Masses:")
        print(self.masses)
        print("Coupling:")
        print(self.coupling)
        print("Time:")
        print(self.time)

    def copy(self):
        # Create a new instance with the same attribute values
        new_molecule = Molecule(
            symbols=self.symbols.copy(),
            coordinates=self.coordinates.copy(),
            momenta=self.momenta.copy() if self.momenta is not None else None,
            scf_energy=self.scf_energy.copy() if self.scf_energy is not None else None,
            forces=self.forces.copy() if self.forces is not None else None,
            amplitudes=self.amplitudes.copy() if self.amplitudes is not None else None,
            timestep=self.timestep,
            multiplicity=self.multiplicity,
            dissociation_flags=self.dissociation_flags,
            elecinfo=self.elecinfo,
            masses=self.masses.copy() if self.masses is not None else None,
            coupling=self.coupling.copy() if self.coupling is not None else None,
            time=self.time.copy() if self.time is not None else None
        )
        
        return new_molecule
    
    def to_dict(self):
        # Convert the Molecule instance to a dictionary
        return {
            'Symbols': self.symbols,
            'Coordinates': self.coordinates.tolist(),
            'Momenta': self.momenta.tolist() if self.momenta is not None else None,
            'SCF Energy': self.scf_energy.tolist(),
            'Forces': self.forces.tolist() if self.forces is not None else None,
            'Amplitudes': [(float(a.real), float(a.imag)) for a in self.amplitudes],
            'Timestep': self.timestep,
            'Multiplicity': self.multiplicity,
            'Dissociation Flags': self.dissociation_flags,
            'Masses': self.masses.tolist() if self.masses is not None else None,
            'Time': self.time.tolist() if self.time is not None else None,
            'Coupling': self.coupling.tolist() if self.coupling is not None else None,
        }

    def to_json(self, filename):
        # Serialize the instance to a JSON file
        with open(filename, 'w') as json_file:
            json.dump(self.to_dict(), json_file, indent=2, cls=NumpyEncoder)

    @classmethod
    def from_dict(cls, data):
        amplitudes = data['Amplitudes']
        if isinstance(amplitudes[0], (float, int)):
            # New format: a list of floats
            pass  # You can keep it as it is
        elif isinstance(amplitudes[0], list):
            # Assume it's a list of real and imaginary parts
            amplitudes = [complex(real, imag) for real, imag in amplitudes]

        return cls(
            symbols=data['Symbols'],
            coordinates=data['Coordinates'],
            momenta=data['Momenta'] if 'Momenta' in data and data['Momenta'] is not None else None,
            scf_energy=data['SCF Energy'] if 'SCF Energy' in data and data['SCF Energy'] is not None else None,
            forces=data['Forces'] if 'Forces' in data and data['Forces'] is not None else None,
            amplitudes=np.array(amplitudes),
            timestep=data['Timestep'],
            multiplicity=data['Multiplicity'],
            dissociation_flags=data['Dissociation Flags'],
            masses=data['Masses'] if 'Masses' in data and data['Masses'] is not None else None,
            time = data['Time'] if 'Time' in data and data['Time'] is not None else None, 
            coupling=data['Coupling'] if 'Coupling' in data and data['Coupling'] is not None else None
        )

    @classmethod
    def from_json(cls, filename):
        # Create a new Molecule instance from a JSON file
        with open(filename, 'r') as json_file:
            data = json.load(json_file)
            return cls.from_dict(data)

def create_empty_molecule(natoms):
    """
    Creates a molecule instance with default attributes.

    Parameters
    ----------
    natoms : int
        Number of atoms in the molecule.
    nst : int
        Number of states for SCF energy (currently irrelevant).
    spin_flip : int
        Indicator for inclusion of Spin-Flip DFT (1=yes/0=NO).

    Returns
    -------
    Molecule
        A Molecule instance initialized with default values.
    """
    symbols = [''] * natoms
    coordinates = np.zeros((natoms, 3), dtype=np.float64)
    scf_energy = np.zeros((gv.num_states), dtype=np.float64)
    momenta = np.zeros((natoms, 3), dtype=np.float64)
    forces = np.zeros((natoms, 3,gv.num_states), dtype=np.float64)
    if gv.spin_flip == 1:
        multiplicity = 5
    elif gv.spin_flip == 0:
        multiplicity = 3
    amplitudes = np.zeros((gv.num_states), dtype=np.complex256)
    amplitudes[0] = 1 + 0j
    masses = np.zeros(natoms, dtype=np.float64)
    return Molecule(symbols, coordinates, momenta, scf_energy, forces, amplitudes, multiplicity=multiplicity, masses=masses)

def initialize_structure():
    """
    Reads molecular structure data from a file and initializes a molecule instance.

    Parameters
    ----------
    nst : int
        Number of SCF states.
    spin_flip : int
        Indicator for inclusion of Spin-Flip DFT (1=yes/0=NO).
    mult : int
        Initial multiplicity setting.

    Functions
    ---------
    
    setup_masses - init.py

    Returns
    -------
    Molecule
        A Molecule instance with attributes initialized based on file data.
    """ 
    file_path = "Geometry"  # Updated file path

    with open(file_path, 'r') as file:
        lines = file.readlines()
    natoms = 0
    for i, line in enumerate(lines):
        if line.strip().lower() == 'momentum':
            natoms = i
            break

    if natoms == 0:
        raise ValueError("The file does not contain a 'momentum' line.")
    
    # Read geometry coordinates
    geometry_lines = lines[0:natoms]
    symbols = [line.split()[0] for line in geometry_lines]
    geometry_data = np.array([list(map(float, line.split()[1:])) for line in geometry_lines], dtype=np.float64)

    # Read momentum coordinates
    momentum_lines = lines[natoms+1:2*natoms+1]
    momentum_data = np.array([list(map(float, line.split())) for line in momentum_lines], dtype=np.float64)

    # Read amplitudes if present
    amplitudes = np.zeros((gv.num_states), dtype=np.complex128)
    amplitudes[gv.start_state-1] = 1 + 0j

    if gv.spin_flip == 1:
        multiplicity = gv.multiplicity+2
    elif gv.spin_flip == 0:
        multiplicity = gv.multiplicity

    forces = np.zeros((natoms, 3, gv.num_states), dtype=np.float64)
    scf_energy = np.zeros((gv.num_states), dtype=np.float64)
    masses = np.zeros(natoms, dtype=np.float64)  # Initialize masses array
    masses = setup_masses(symbols)
    coupling = np.zeros((natoms,3))
    time = np.zeros((4))
    # Create Molecule object with default amplitudes
    molecule = Molecule(symbols, geometry_data, momenta=momentum_data, scf_energy=scf_energy, forces=forces, multiplicity=multiplicity, amplitudes=amplitudes, masses=masses,coupling=coupling,time=time)

    return molecule

def create_molecule(natoms = None):
    """
    Wrapper function for molecule creation with optional structure reading.

    This function utilizes the `create_empty_molecule` function to create 
    a new Molecule instance with default attributes when `reps` is None. 
    If `reps` is provided, it calls `initialize_structure` to read data 
    from the "Geometry" file to initialize the molecule attributes.


    Parameters
    ----------
    reps : any
        Input to determine whether to initialize structure from a file.
    natoms : int
        Number of atoms in the molecule.
    nst : int
        Number of SCF states.
    spin_flip : int
        Spin flip indicator.

    Functions
    ---------

    create_empty_molecule - init.py

    initialise_structure - init.py    
    
    Returns
    -------
    Molecule
        A Molecule instance created with either default or file-based attributes.
    """
    if natoms is not None:
        # If reps is None, create an empty molecule
        return create_empty_molecule(natoms)
    else:
        # Otherwise, create an empty molecule
        return initialize_structure()

def setup_masses(symbols):
    """
    Returns atomic masses for a given set of symbols and converts from atomic mass to atomic mass units (amu).

    Parameters
    ----------
    symbols : list[str]
        List of chemical symbols representing the molecule.

    Returns
    -------
    np.ndarray
        Array of atomic masses for the specified symbols in amu.
    """
    ATOMIC_MASSES = {
        'H': 1, 'He': 4, 'Li': 7, 'Be': 9, 'B': 11, 'C': 12, 'N': 14, 'O': 16, 'F': 19, 'Ne': 20,
        'Na': 23, 'Mg': 24, 'Al': 27, 'Si': 28, 'P': 31, 'S': 32, 'Cl': 35, 'Ar': 40, 'K': 39, 'Ca': 40,
        'Sc': 45, 'Ti': 48, 'V': 51, 'Cr': 52, 'Mn': 55, 'Fe': 56, 'Co': 59, 'Ni': 59, 'Cu': 64, 'Zn': 65,
        'Ga': 70, 'Ge': 73, 'As': 75, 'Se': 79, 'Br': 80, 'Kr': 84, 'Rb': 85, 'Sr': 88, 'Y': 89, 'Zr': 91,
        'Nb': 93, 'Mo': 96, 'Tc': 98, 'Ru': 101, 'Rh': 103, 'Pd': 106, 'Ag': 108, 'Cd': 112, 'In': 115,
        'Sn': 119, 'Sb': 122, 'Te': 127, 'I': 127, 'Xe': 131, 'Cs': 133, 'Ba': 137, 'La': 138, 'Ce': 140,
        'Pr': 141, 'Nd': 144, 'Pm': 145, 'Sm': 150, 'Eu': 152, 'Gd': 157, 'Tb': 159, 'Dy': 162, 'Ho': 165,
        'Er': 167, 'Tm': 169, 'Yb': 173, 'Lu': 175, 'Hf': 178, 'Ta': 181, 'W': 184, 'Re': 186, 'Os': 192,
        'Ir': 193, 'Pt': 195, 'Au': 197, 'Hg': 201, 'Tl': 204, 'Pb': 207, 'Bi': 209, 'Th': 232, 'U': 238
    }

    
    # Fetch the masses from the dictionary
    masses = np.array([ATOMIC_MASSES.get(symbol, 0.0) for symbol in symbols], dtype=np.float64)
    
    # Convert to mau units (1 amu = 1822.887 mau)
    masses_mau = masses * 1822.887
    
    return masses_mau

def check_timestep_file(filename):
    """Check if the simulation can be restarted based on the existing XYZ file."""
    timesteps = []
    if os.path.exists(filename):
        with open(filename, "r") as xyz_file:
            lines = xyz_file.readlines()
        for line in reversed(lines):
            if "Timestep:" in line:
                line = line.split()
                timesteps.append(int(line[1]))
    return timesteps

def check_timestep_file(filename):
    """Check if the simulation can be restarted based on the existing XYZ file."""
    timesteps = []
    if os.path.exists(filename):
        with open(filename, "r") as xyz_file:
            lines = xyz_file.readlines()
        for line in reversed(lines):
            if "Timestep:" in line:
                line = line.split()
                timesteps.append(int(line[1]))
    return timesteps

def read_xyz_file(filename, target_timestep):
    """Read XYZ file and return symbols and coordinates for the target timestep."""
    symbols = []
    coordinates = []
    reading = False
    if os.path.exists(filename):
        with open(filename, "r") as xyz_file:
            lines = xyz_file.readlines()
        for line in lines:
            if "Timestep:" in line:
                current_timestep = int(line.split()[1])
                if current_timestep == target_timestep:
                    reading = True
                else:
                    reading = False
            elif reading:
                if line.strip() == "----------------------------------------":
                    break
                parts = line.split()
                if len(parts) >= 5:
                    symbols.append(parts[1])
                    coordinates.append([float(parts[2]), float(parts[3]), float(parts[4])])
    return symbols, np.array(coordinates)

def read_momenta_file(filename, target_timestep):
    """Read Momenta file and return momenta for the target timestep."""
    momenta = []
    reading = False
    if os.path.exists(filename):
        with open(filename, "r") as momenta_file:
            lines = momenta_file.readlines()
        for line in lines:
            if "Timestep:" in line:
                current_timestep = int(line.split()[1])
                if current_timestep == target_timestep:
                    reading = True
                else:
                    reading = False
            elif reading:
                if line.strip() == "----------------------------------------":
                    break
                parts = line.split()
                if len(parts) == 3:
                    momenta.append([float(parts[0]), float(parts[1]), float(parts[2])])
    return np.array(momenta)

import os
import numpy as np

def read_forces_file(filename, target_timestep):
    """Read Forces file and return forces, amplitudes, multiplicity, SCF energy, and optional coupling forces for the target timestep."""

    forces_dict = {}  # Store forces by atom and state
    coupling_dict = {}  # Store coupling forces separately
    amplitudes = None
    multiplicity = None
    scf_energy = None
    dissociation_flags = None
    reading = False
    current_atom = None
    current_state = None
    num_states = 1  # Start at 1, assuming at least one state exists

    if os.path.exists(filename):
        with open(filename, "r") as forces_file:
            lines = forces_file.readlines()

        for line in lines:
            line = line.strip()

            # Detect timestep and check if we should start reading
            if line.startswith("Timestep:"):
                try:
                    current_timestep = int(line.split()[1])
                    reading = current_timestep == target_timestep
                except (ValueError, IndexError):
                    reading = False
                if reading:
                    forces_dict.clear()
                    coupling_dict.clear()
                    num_states = 1  # Reset for new timestep

            elif reading:
                if line == "----------------------------------------":
                    break  # Stop reading at the end of the timestep block

                # Parse amplitudes
                elif line.startswith("Amplitudes:"):
                    amplitude_str = line.replace("Amplitudes:", "").strip().strip("[]")
                    amplitudes = np.array([complex(x) for x in amplitude_str.split()])

                # Parse multiplicity
                elif line.startswith("Multiplicity:"):
                    try:
                        multiplicity = int(line.split()[1])
                    except (ValueError, IndexError):
                        multiplicity = None

                # Parse SCF energy
                elif line.startswith("SCF Energy:"):
                    energy_str = line.replace("SCF Energy:", "").strip().strip("[]")
                    try:
                        scf_energy = np.array([float(x) for x in energy_str.split()])
                    except ValueError:
                        scf_energy = None

                # Parse dissociation flags
                elif line.startswith("Dissociation Flags:"):
                    flags_str = line.replace("Dissociation Flags:", "").strip().strip("[]")
                    dissociation_flags = [flag.strip().strip("'") for flag in flags_str.split(",")]


                # Parse forces per atom and state
                elif line.startswith("Atom"):
                    try:
                        current_atom = int(line.split()[1].strip(":"))  # Extract atom number
                        forces_dict[current_atom] = {}  # Initialize forces for the atom
                    except (ValueError, IndexError):
                        current_atom = None

                elif line.startswith("State") and current_atom is not None:
                    try:
                        parts = line.split()
                        current_state = int(parts[1].strip(":"))  # Extract state number
                        num_states = max(num_states, current_state)  # Ensure we track max state index
                        forces_dict[current_atom][current_state] = np.array(
                            [float(parts[2]), float(parts[3]), float(parts[4])]
                        )

                    except (ValueError, IndexError):
                        pass  # Skip if there is an issue parsing the line

                # Parse coupling forces
                elif line.startswith("Coupling:") and current_atom is not None:
                    try:
                        parts = line.split()
                        coupling_dict[current_atom] = np.array(
                            [float(parts[1]), float(parts[2]), float(parts[3])]
                        )
                    except (ValueError, IndexError):
                        pass  # Skip if there is an issue parsing the coupling forces

    # Convert forces into a structured NumPy array
    if forces_dict:
        natoms = max(forces_dict.keys())  # Total number of atoms
        num_states = max(max(states.keys()) for states in forces_dict.values())  # Ensure num_states is correctly determined
        num_states = max(num_states, 1)  # Ensure at least one state


        forces = np.zeros((natoms, 3, num_states))  # Allocate (natoms, 3, num_states)



        for atom, states in forces_dict.items():
            for state, force in states.items():   
                forces[atom - 1, :, state - 1] = force


    else:
        forces = None

    # Convert coupling forces into a structured NumPy array
    if coupling_dict:
        natoms = max(coupling_dict.keys()) if forces_dict else 0
        coupling_forces = np.zeros((natoms, 3))
        for atom, coupling in coupling_dict.items():
            coupling_forces[atom - 1, :] = coupling  # Adjust for zero-based indexing
    else:
        coupling_forces = None  # If no coupling terms are present

    return forces, amplitudes, multiplicity, scf_energy, dissociation_flags, coupling_forces


def read_times_file(filename, target_timestep):
    """Read times from the time.all file for the target timestep."""
    times = None
    reading = False
    if os.path.exists(filename):
        with open(filename, "r") as times_file:
            lines = times_file.readlines()
        for line in lines:
            if reading == False and "Timestep:" in line:
                current_timestep = int(line.split()[1])
                reading = current_timestep == target_timestep
            elif reading:
                if line.startswith("Time of this step:"):
                    time_of_step = float(line.split(":")[1].strip().rstrip("s"))
                elif line.startswith("Time of QChem this step:"):
                    time_qchem = float(line.split(":")[1].strip().rstrip("s"))
                elif line.startswith("Total time of propagation:"):
                    total_time = float(line.split(":")[1].split()[0].strip().rstrip("s"))
                elif line.startswith("Total time of QChem:"):
                    total_time_QChem = float(line.split(":")[1].split()[0].strip().rstrip("s"))
                                # Now break only if at least one time value was found
                elif line.strip() == "----------------------------------------":
                    if time_of_step is not None:
                        break  
        times = np.array([time_qchem, time_of_step, total_time_QChem, total_time])
    return times

def find_timestep():
    """Check and return the congruent timestep from various files."""
    try:
        with open("output/molecule.json") as f:
            molecule_json = json.load(f)
        
        xyz_timesteps = check_timestep_file('output/xyz.all')
        momenta_timesteps = check_timestep_file('output/momenta.all')
        forces_timesteps = check_timestep_file('output/forces.all')
        json_timestep = molecule_json.get('Timestep', 0)
        for timestep in xyz_timesteps:
            if timestep in momenta_timesteps and timestep in forces_timesteps:
                if timestep == json_timestep:
                    start_step = timestep / gv.timestep
                    molecule1 = Molecule.from_json('output/molecule.json')
                    print('coordinates', molecule1.coordinates)
                    num_atoms = len(molecule1.symbols)
                    molecule2 = create_empty_molecule(num_atoms)
                    print('Found congruent timestep: ', timestep, 'so start_step is:', start_step)
                    return molecule1, molecule2, start_step
                else:
                    print('JSON file disagrees. Reading molecule from XYZ, Momenta, and Forces files.')
                    symbols, coordinates = read_xyz_file('output/xyz.all', timestep)
                    print("coordinates: ", coordinates)
                    momenta = read_momenta_file('output/momenta.all', timestep)
                    forces, amplitudes, multiplicity, scf_energy, dissociation_flags, coupling = read_forces_file('output/forces.all', timestep)
                    time = read_times_file('output/time.all', timestep)
                    masses = setup_masses(symbols)
                    molecule1 = Molecule(symbols, coordinates, momenta, scf_energy, forces, amplitudes, timestep, multiplicity, dissociation_flags=dissociation_flags, time = time, coupling = coupling, masses = masses)
                    num_atoms = len(molecule1.symbols)
                    molecule2 = create_empty_molecule(num_atoms)
                    start_step = timestep / gv.timestep
                    return molecule1, molecule2, start_step 

        print('No congruent timestep found. Setting start_step to 0.')
        return 0

    except FileNotFoundError as e:
        print(f"Error: {e}")
        return 0

def initialize_simulation():
    """Check if restart is possible and determine the start step."""
    # Check if the result files exist
    result_files_exist = all(os.path.exists(f'output/{file}') for file in ['dissociation.out', 'fragments.out'])

    if result_files_exist:
        print('Result files exist. Stopping the simulation.')
        exit(0)

    # Check if the necessary files exist for restart
    files_exist = all(os.path.exists(f'output/{file}.all') for file in ['xyz', 'momenta', 'forces']) and os.path.exists('output/molecule.json')

    if files_exist:
        # Find the latest congruent timestep
        molecule1, molecule2, start_step = find_timestep()
    else:
        # If files don't exist, set start_step to 0
        start_step = 0
        molecule1 = initialize_structure()
        num_atoms = len(molecule1.symbols)
        molecule2 = create_empty_molecule(num_atoms)
        molecule1 = elec.run_elec_structure(molecule1, Guess=False)
        

    return molecule1, molecule2, start_step
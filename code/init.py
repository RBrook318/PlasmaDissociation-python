import numpy as np
import json
np.set_printoptions(precision=30)

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.float128):
            return float(obj)
        return super(NumpyEncoder, self).default(obj)

class Molecule:
    def __init__(self, symbols, coordinates, momenta=None, scf_energy=None, forces=None, amplitudes=None, timestep=0, multiplicity=5, dissociation_flags=None, elecinfo=None, masses=None):
        self.symbols = symbols
        self.coordinates = np.array(coordinates, dtype=np.float128)
        self.momenta = np.array(momenta, dtype=np.float128) if momenta is not None else None
        self.scf_energy = np.array(scf_energy, dtype=np.float128)
        self.forces = np.array(forces, dtype=np.float128) if forces is not None else None
        self.amplitudes = np.array(amplitudes, dtype=np.complex256) if amplitudes is not None else np.array([1.0 + 0.0j], dtype=np.complex256)
        self.timestep = timestep
        self.multiplicity = multiplicity
        self.dissociation_flags = dissociation_flags or ["NO"] * len(symbols)
        self.elecinfo = elecinfo if elecinfo is not None else None
        self.masses = np.array(masses, dtype=np.float128) if masses is not None else np.zeros(len(symbols), dtype=np.float128)

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
        self.masses = np.array(new_masses, dtype=np.float128)

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
            masses=self.masses.copy() if self.masses is not None else None
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
            'Masses': self.masses.tolist() if self.masses is not None else None
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
            masses=data['Masses'] if 'Masses' in data and data['Masses'] is not None else None
        )

    @classmethod
    def from_json(cls, filename):
        # Create a new Molecule instance from a JSON file
        with open(filename, 'r') as json_file:
            data = json.load(json_file)
            return cls.from_dict(data)


def create_empty_molecule(natoms, nst, spin_flip):
    symbols = [''] * natoms
    coordinates = np.zeros((natoms, 3), dtype=np.float128)
    scf_energy = np.zeros((nst), dtype=np.float128)
    momenta = np.zeros((natoms, 3), dtype=np.float128)
    forces = np.zeros((natoms, 3), dtype=np.float128)
    if spin_flip == 1:
        multiplicity = 5
    elif spin_flip == 0:
        multiplicity = 3
    amplitudes = np.zeros((nst), dtype=np.complex256)
    amplitudes[0] = 1 + 0j
    masses = np.zeros(natoms, dtype=np.float128)
    return Molecule(symbols, coordinates, momenta, scf_energy, forces, amplitudes, multiplicity=multiplicity, masses=masses)

def initialize_structure(nst, spin_flip):
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
    if len(lines) > 2*natoms+1:
        amplitudes_lines = lines[2*natoms+2:]  # Assuming amplitudes are two lines below "momentum"
        amplitudes = list(map(float, amplitudes_lines))
    else: 
        amplitudes = [1.0, 0.0]

    if spin_flip == 1:
        multiplicity = 5
    elif spin_flip == 0:
        multiplicity = 3

    forces = np.zeros((natoms, 3), dtype=np.float128)
    scf_energy = np.zeros((nst), dtype=np.float128)
    masses = np.zeros(natoms, dtype=np.float128)  # Initialize masses array
    masses = setup_masses(symbols)
    # Create Molecule object with default amplitudes
    molecule = Molecule(symbols, geometry_data, momenta=momentum_data, scf_energy=scf_energy, forces=forces, multiplicity=multiplicity, amplitudes=amplitudes, masses=masses)

    return molecule

def create_molecule(reps, natoms, nst, spin_flip):
    if reps is None:
        # If reps is None, create an empty molecule
        return create_empty_molecule(natoms, nst, spin_flip)
    else:
        # Otherwise, create an empty molecule
        return initialize_structure(nst, spin_flip)

def setup_masses(symbols):
    
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
    masses = np.array([ATOMIC_MASSES.get(symbol, 0.0) for symbol in symbols], dtype=np.float128)
    
    # Convert to mau units (1 amu = 1822.887 mau)
    masses_mau = masses * 1822.887
    
    return masses_mau
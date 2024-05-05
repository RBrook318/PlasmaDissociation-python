import numpy as np

def calculate_distance(coord1, coord2):
    return np.linalg.norm(coord1 - coord2)

def print_distances(coords_list, atom_index):
    for i, coord in enumerate(coords_list):
        if i != atom_index:
            distance = calculate_distance(coords_list[atom_index], coord)
            print(f"Distance from atom {atom_index} to atom {i}: {distance}")

# Example usage
atoms_data = []
with open("xyz.txt", "r") as file:
    for line in file:
        parts = line.strip().split()
        coords = np.array([float(parts[1]), float(parts[2]), float(parts[3])])
        atoms_data.append(coords)

# Print distances from atom 0 to all other atoms
atom_index = 13
print_distances(atoms_data, atom_index)
import numpy as np

def parse_xyz_all(file_path):
    timesteps_data = []
    atoms_data = []

    with open(file_path, 'r') as file:
        timestep_atoms = []
        for line in file:
            line = line.strip()
            if line.startswith("Timestep:"):
                if timestep_atoms:
                    timesteps_data.append(timestep_atoms)
                    timestep_atoms = []
            elif line.startswith("-"):
                continue
            elif line:
                parts = line.split()
                coords = np.array([float(parts[2]), float(parts[3]), float(parts[4])])
                element = parts[1]
                atom = {
                    'element': element,
                    'coords': coords
                }
                atoms_data.append(atom)
                timestep_atoms.append(atom)

        # Append the last timestep
        if timestep_atoms:
            timesteps_data.append(timestep_atoms)

    return timesteps_data, atoms_data
timesteps_data, atoms_data = parse_xyz_all("xyz.all")
def calculate_distances(atoms):
    num_atoms = len(atoms)
    distances = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = np.linalg.norm(atoms[i]['coords'] - atoms[j]['coords'])
            distances[i][j] = distances[j][i] = dist

    return distances

def find_fragments(atoms, distances, threshold=5):
    fragments = []
    visited = set()

    def dfs(atom_index, fragment):
        visited.add(atom_index)
        fragment.append(atom_index)
        for neighbor_index, distance in enumerate(distances[atom_index]):
            if neighbor_index not in visited and distance <= threshold:
                dfs(neighbor_index, fragment)

    for i in range(len(atoms)):
        if i not in visited:
            fragment = []
            dfs(i, fragment)
            fragments.append(fragment)

    return fragments

def is_fully_separated(fragment, all_fragments, distances, threshold=5):
    for other_fragment in all_fragments:
        if fragment != other_fragment:
            for atom_index in fragment:
                for other_atom_index in other_fragment:
                    if distances[atom_index][other_atom_index] <= threshold:
                        return False
    return True

threshold = 30
old_length = 0

with open("isolated_fragments_output.txt", "w") as output_file:
    for timestep, atoms in enumerate(timesteps_data):
        distances = calculate_distances(atoms)
        fragments = find_fragments(atoms, distances, threshold)
        fully_separated_fragments = []
        
        for fragment in fragments:
            if is_fully_separated(fragment, fragments, distances, threshold):
                fully_separated_fragments.append(fragment)
        
        length = len(fully_separated_fragments)
        if length > old_length:
            output_file.write(f"Timestep: {timestep}\n")
            output_file.write("Isolated Fragments:\n")
            output_file.write(f"{fully_separated_fragments}\n")

            # Print distances between atoms from different fragments
            output_file.write("Distances between atoms from different fragments:\n")
            for i, fragment1 in enumerate(fully_separated_fragments):
                for j, fragment2 in enumerate(fully_separated_fragments):
                    if i != j:
                        for atom_index1 in fragment1:
                            for atom_index2 in fragment2:
                                distance = distances[atom_index1][atom_index2]
                                output_file.write(f"Distance from {atom_index1+1} inf{i+1} "
                                                   f"and  {atom_index2+1} inf {j+1}: {distance}\n")
            output_file.write("\n")
            
            old_length = length

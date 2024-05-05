# code for splitting the molecule into different parts. To be called from main in order to check for isolated groups. 
# This should be done by checking for isolate groups within the molecule, then it should separate those into molecule objects 
# to be able to separate their qchem runs and then recombined for the output
import init
import numpy as np

def full_isolates(molecule_array, nisolates):

    if nisolates==1:
        distances = calculate_distances(molecule_array[0])
        fragments = find_fragments(molecule_array[0], distances, 30)
        fully_separated_fragments = []
        
        for fragment in fragments:
            if is_fully_separated(fragment, fragments, distances, 30):
                modified_fragment = [atom_index + 1 for atom_index in fragment]
                fully_separated_fragments.append(modified_fragment)
        print(fully_separated_fragments)

        if len(fully_separated_fragments)>1: 
            for i in range(1,len(fully_separated_fragments)):
                sub_molecule,parent_molecule = init.sub_molecule(molecule_array[0], fully_separated_fragments[i])
                print('parent molecule: ')
                parent_molecule.print_info()
                print('sub molecule: ')
                sub_molecule.print_info()
                molecule_array.append(parent_molecule)
                molecule_array.append(sub_molecule)
    else: 
        for i in range(1,nisolates+1): # go through each sub_molecule in the array 
            distances = calculate_distances(molecule_array[i])
            fragments = find_fragments(molecule_array[i], distances, 30)
            fully_separated_fragments = [] 
            
            for fragment in fragments:
                if is_fully_separated(fragment, fragments, distances, 30):
                    modified_fragment = [atom_index + 1 for atom_index in fragment]
                    fully_separated_fragments.append(modified_fragment)
            print(fully_separated_fragments)

            if len(fully_separated_fragments)>1: # go through each fragment found
                for j in range(1,len(fully_separated_fragments)+1):
                    sub_molecule,parent_molecule = init.sub_molecule(molecule_array) # The parent molecule will be the moelcule where the fragment was found. I need to make some thing to separate the different sub molecules from the larger fragment it was found in. 
                    # Then the parent molecule replaces where it used to be in the array of molecules and the new generated one will go at the end. As there is no interaction between them, there is no need to keep a record of their connection. Great news.
                    # Just now need to code the idea of sub molecules, the recombination and then the multiple qchem submits and retrievals. 
                    molecule_array.append(sub_molecule)
    
    nisolates=len(molecule_array)
    
    return molecule_array, nisolates


def calculate_distances(molecule):
    num_atoms = len(molecule.symbols)
    distances = np.zeros((num_atoms, num_atoms))

    for i in range(num_atoms):
        for j in range(i + 1, num_atoms):
            dist = np.linalg.norm(molecule.coordinates[i]- molecule.coordinates[j])
            distances[i][j] = distances[j][i] = dist

    return distances

def find_fragments(molecule, distances, threshold=5):
    fragments = []
    visited = set()

    def dfs(atom_index, fragment):
        visited.add(atom_index)
        fragment.append(atom_index)
        for neighbor_index, distance in enumerate(distances[atom_index]):
            if neighbor_index not in visited and distance <= threshold:
                dfs(neighbor_index, fragment)

    for i in range(len(molecule.symbols)):
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

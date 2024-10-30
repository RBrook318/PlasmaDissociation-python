# Result file for the pure python plasma code. Should use a provided bondingarray file to know which bonds to investigate. Then it can look through xyz.all and 
# then make a dissociation.out as well as the C-C.out order files both in the output folder and in the total folder of the a whole set of trajectories. 
import os 
import json
from init import Molecule
import networkx as nx
import shutil
import specifics as spc
from collections import defaultdict
import graphs as grph

def process_results():
    # 0. Increment count file
    track_completed_trajectory()
    # 1.read in the bonding array
    bondarr = read_bondarr()
    # # 2. detect dissociation
    detect_dissociation(bondarr)
    # 3. Find fragments.
    fragments()
    # 4. compile dissociation data
    compile_results()
    # 5. Compile fragment data
    combine_fragments()
    # 6. Specifics
    spc.specifics('../results/bondarr.txt')
    #  7. graphs
    write_graphs_json()
    grph.create_graphs()

def read_bondarr():
    """
    Load the bonding array from a specified file.

    Parameters:
    ----------
    file_path : str
        Path to the bonding array file containing bond information.

    Returns:
    -------
    bonds : list
        A list of bond entries, each representing a bond event in the array.
    """
    bondarr = {}
    with open('../results/bondarr.txt', 'r') as file:
        for line in file:
            line = line.strip()
            if line:  # Skip empty lines
                atoms, bond_type = line.split(':')
                atom1, atom2 = atoms.split('-')
                atom1 = int(atom1)
                atom2 = int(atom2)
                bondarr[(atom1, atom2)] = bond_type
    
    return bondarr

def detect_dissociation(bondarr):

    """
    Detects bond dissociations in molecular dynamics data based on atomic distances.

    This function reads atomic positions over time from the `xyz.all` file, iterating 
    through each timestep and calculating the distances between specified bonded pairs 
    of atoms. If a bond's distance exceeds a defined threshold (5.0 units), it logs the 
    dissociation event along with the bond type in the `dissociation.out` file. 
    Each bond dissociation is recorded only once per simulation run to prevent duplicate entries.

    Parameters:
    -----------
    bondarr : dict
        A dictionary where each key is a tuple of atom indices (i, j) representing a bonded pair, 
        and the value is a string describing the bond type.

    Outputs:
    --------
    Writes bond dissociation events to `dissociation.out` with details on the timestep, 
    broken bond, and bond type.

    Notes:
    ------
    - The dissociation threshold is set at 5.0 units by default.
    - Only dissociations not previously recorded are logged to avoid duplicates.

    Example Entry in `dissociation.out`:
    ------------------------------------
    "Dissociation detected at timestep 100.0, Broken bond: 1-2:C-C"
    """
    # Read the input data from a file
    with open('output/xyz.all', 'r') as f:
        lines = f.readlines()

    # Initialize variables to store the current timestep, atom data, and dissociation flag
    timestep = None
    atoms = {}
    dissociated_bonds = set()  # Keep track of dissociated bond pairs

    # Open the output file for writing
    with open('output/dissociation.out', 'w') as output:
        # Iterate through the lines and process the data
        for line in lines:
            parts = line.split()
            if len(parts) == 2 and parts[1].isdigit():
                if atoms:
                    for bonded_pair, bond_type in bondarr.items(): 
                        if bonded_pair not in dissociated_bonds: 
                            i, j = bonded_pair
                            atom1 = atoms[i]        
                            atom2 = atoms[j]
                            distance = ((atom1[1] - atom2[1])**2 + (atom1[2] - atom2[2])**2 + (atom1[3] - atom2[3])**2)**0.5
                            if distance > 5.0:  # Adjust this threshold as needed
                                broken_bond_str = f"{i}-{j}:{bond_type}"  # Include bond type in the output
                                if broken_bond_str not in dissociated_bonds:
                                    output.write(f"Dissociation detected at timestep {timestep}, Broken bond: {broken_bond_str}\n")
                                    dissociated_bonds.add(bonded_pair)
                                    
                    atoms = {}
                timestep = float(parts[1])
            elif len(parts) == 5 and parts[0].isdigit():
                atom_num = int(parts[0])
                atom_info = [parts[1]] + [float(x) for x in parts[2:]]
                atoms[atom_num] = atom_info

def compile_results():
    """
    Compiles dissociation results from the `dissociation.out` file and logs them into 
    a summary file along with bond-specific files.

    This function reads dissociation events from the `dissociation.out` file, checks for 
    any detected dissociations, and appends the results to a collated summary file 
    (`collated_diss.txt`). It also generates bond-type-specific output files for each 
    bond type found in the dissociation events and a general file for all bonds.

    The results are structured as follows:
    - If dissociations are detected, the details are appended to the summary file with 
      the current run number.
    - If no dissociations are found, an entry indicating this is logged in the summary file.
    - For each dissociation event, the timestep is logged in both bond-type-specific files 
      and an "all bonds" file.

    Outputs:
    --------
    - Appends to `../results/collated_diss.txt` to summarize dissociation events.
    - Creates or appends to bond-type-specific output files in `../results/bonds/` 
      for each detected bond type and an "allbonds.out" file.

    Notes:
    ------
    - The run number is extracted from the current working directory.
    - Each output file will be appended to, preserving previous entries.

    Example Entry in `collated_diss.txt`:
    ---------------------------------------
    "--- run-1 ---\n"
    "Dissociation detected at timestep 100.0, Broken bond: 1-2:C-C\n"
    "----------------------------------------\n"

    Example Entry in bond-type-specific file:
    -------------------------------------------
    "100\n"  # Timestep indicating when the bond was broken
    """
    
    with open('output/dissociation.out', "r") as f:
        lines = f.readlines()
    EXDIR= os.getcwd()
    print(EXDIR)
    rep_number = EXDIR.split('/')[-1].split('-')[-1]
    if lines:  # Only proceed if lines is not empty
        with open('../results/collated_diss.txt', 'a') as output_file:
            output_file.write(f'--- run-{rep_number} ---\n')
            output_file.write(''.join(lines))  # Join the lines into a single string
            output_file.write('----------------------------------------\n')
    else: 
        with open('../results/collated_diss.txt', 'a') as output_file:
            output_file.write(f'--- run-{rep_number} ---\n')
            output_file.write('No dissociation found \n')  # Join the lines into a single string
            output_file.write('----------------------------------------\n')

    for line in lines:
        if line.startswith("Dissociation detected"):
            parts = line.split(", ")
            timestep = int(float(parts[0].split(" ")[-1]))
            bond_info = parts[1].split(":")
            bond_type = bond_info[2].strip() # Remove leading and trailing whitespace
            bond_number = bond_info[1].strip()
            
            # Write to bond-type-specific output file
            bond_type_file = os.path.join("../results/bonds/", f"{bond_type}.out")
            with open(bond_type_file, "a") as f_out:
                f_out.write(f"{timestep}\n")
                
            # Write to bond-type-specific output file
            bond_type_file = os.path.join("../results/bonds/", "allbonds.out")
            with open(bond_type_file, "a") as f_out:
                f_out.write(f"{timestep}\n")

            # Write to old output file format and order by timestep
            bond_number_file = os.path.join("../results/bonds/", f"{bond_number}.out")
            with open(bond_number_file, "a") as f_out:
                f_out.write(f"{timestep}\n")


def fragments():
    """
    Analyzes the molecular structure from a JSON file to identify and categorize fragments 
    based on atomic connectivity and distances.

    This function reads atomic coordinates and symbols from a molecule JSON file, calculates 
    distances between atoms, and constructs a graph to identify connected fragments. It then 
    determines the molecular formula for each fragment and logs the results into an output file.

    The process involves the following steps:
    1. Load the molecule data from 'output/molecule.json'.
    2. Calculate pairwise distances between atoms and create a connectivity graph.
    3. Identify connected components (fragments) in the graph.
    4. For each connected fragment, compute the molecular formula, organized by 
       Carbon (C), Fluorine (F), Hydrogen (H), and Oxygen (O) counts. (SHOULD ADD Si SOON)
    5. Handle isolated atoms as individual fragments and compute their formulas.
    6. Write the fragment formulas and their counts to 'output/fragments.out'.

    Outputs:
    --------
    - Writes to 'output/fragments.out', detailing each unique fragment formula and its count.

    Example Output in `fragments.out`:
    -----------------------------------
    Fragment Formulas:
    C3H8: 2
    C2H4: 1
    O: 1
    F2: 1

    Notes:
    ------
    - Distances less than 5.0 units are considered to indicate a bond between atoms.
    - The molecular formula is constructed in the order of C, F, H, and O.
    """
    molecule= Molecule.from_json('output/molecule.json')
    coordinates = molecule.coordinates
    fragment_formulas = {} 
    atom_coordinates = [(i, molecule.symbols[i], molecule.coordinates[i, 0], molecule.coordinates[i, 1], molecule.coordinates[i, 2]) for i in range(len(molecule.symbols))]
    atom_distances = {}
    for i in range(len(atom_coordinates)):
        for j in range(i + 1, len(atom_coordinates)):
            dist = distance(atom_coordinates[i], atom_coordinates[j])
            atom_distances[(i, j)] = dist

    G = nx.Graph()
    for (i, j), dist in atom_distances.items():
        if dist < 5:
            G.add_edge(i, j)

    # Initialize a set to keep track of used atoms
    used_atoms = set()
    connected_fragments = list(nx.connected_components(G))
    fragments = []

    for idx, fragment in enumerate(connected_fragments):
        fragment_atoms = []
        for atom_idx in fragment:
            atom_number, atom_name, x, y, z = atom_coordinates[atom_idx]
            fragment_atoms.append((atom_number, atom_name, x, y, z))
            used_atoms.add(atom_idx)
        fragments.append(fragment_atoms)

        # Calculate the molecular formula for the fragment
        formula = {}
        for atom in fragment_atoms:
            atom_name = atom[1]
            if atom_name in formula:
                formula[atom_name] += 1
            else:
                formula[atom_name] = 1

        # Organize formula by CFHO (include 'O' in the order)
        cfho_formula = ""
        if 'C' in formula:
            cfho_formula += f"C{formula['C']}"
        if 'F' in formula:
            cfho_formula += f"F{formula['F']}"
        if 'H' in formula:
            cfho_formula += f"H{formula['H']}"
        if 'O' in formula:
            cfho_formula += f"O{formula['O']}"

        if cfho_formula in fragment_formulas:
            fragment_formulas[cfho_formula] += 1
        else:
            fragment_formulas[cfho_formula] = 1

    # Process isolated atoms as separate fragments
    isolated_atoms = set(range(len(atom_coordinates))) - used_atoms
    for atom_idx in isolated_atoms:
        atom_number, atom_name, x, y, z = atom_coordinates[atom_idx]
        fragments.append([(atom_number, atom_name, x, y, z)])

        # Calculate the molecular formula for the isolated atom
        formula = {atom_name: 1}

        # Organize formula by CFH
        cfho_formula = ""
        if 'C' in formula:
            cfho_formula += f"C{formula['C']}"
        if 'F' in formula:
            cfho_formula += f"F{formula['F']}"
        if 'H' in formula:
            cfho_formula += f"H{formula['H']}"
        if 'O' in formula:
            cfho_formula += f"O{formula['O']}"

        if cfho_formula in fragment_formulas:
            fragment_formulas[cfho_formula] += 1
        else:
            fragment_formulas[cfho_formula] = 1

    with open("output/fragments.out", "w") as count_file:
        count_file.write("Fragment Formulas:\n")
        for formula, count in fragment_formulas.items():
            count_file.write(f"{formula}: {count}\n")


def distance(point1, point2):
    x1, y1, z1 = point1[2:]
    x2, y2, z2 = point2[2:]
    return ((x1 - x2) ** 2 + (y1 - y2) ** 2 + (z1 - z2) ** 2) ** 0.5

def combine_fragments():
    if os.path.exists("../results/fragments.out"):
       combine_fragment_counts("output/fragments.out", "../results/fragments.out", "../results/fragments.out")
    else:
        shutil.copy2("output/fragments.out","../results")

def combine_fragment_counts(file1, file2, output_file):
    """
    Combines fragment counts from two files and writes the results to an output file.

    This function reads fragment counts from two specified input files, combines the counts 
    for each unique fragment formula, and writes the combined counts to a specified output file.

    Parameters:
    -----------
    file1 : str
        Path to the first input file containing fragment counts.

    file2 : str
        Path to the second input file containing fragment counts.

    output_file : str
        Path to the output file where combined fragment counts will be written.

    Notes:
    ------
    - Assumes that each input file has a consistent format with fragment formulas and counts 
      listed in a specific manner.
    - Uses helper functions to read counts, combine them, and write the output.
    """
    # Read fragment counts from the first file
    fragment_formulas1 = read_fragment_counts(file1)

    # Read fragment counts from the second file
    fragment_formulas2 = read_fragment_counts(file2)

    # Combine fragment counts from both files
    combined_fragment_formulas = combine_counts(fragment_formulas1, fragment_formulas2)

    # Write the combined fragment counts to the output file
    write_fragment_counts(output_file, combined_fragment_formulas)

def read_fragment_counts(file_path):
    """
    Reads fragment counts from a specified file.

    This function parses the input file to extract fragment formulas and their associated 
    counts, returning them as a dictionary.

    Parameters:
    -----------
    file_path : str
        Path to the input file containing fragment counts.

    Returns:
    --------
    dict
        A dictionary where keys are fragment formulas (str) and values are their counts (int).

    Notes:
    ------
    - The function skips the first line of the file, assuming it contains a header or title.
    - Each subsequent line is expected to be in the format 'formula: count'.
    """
    fragment_formulas = {}
    with open(file_path, 'r') as file:
        lines = file.readlines()

        # Skip the first line (title or header)
        lines = lines[1:]

        for line in lines:
            parts = line.strip().split(':')
            if len(parts) == 2:
                formula, count = parts[0], int(parts[1])
                fragment_formulas[formula] = count
    return fragment_formulas

def combine_counts(counts1, counts2):
    """
    Combines two dictionaries of fragment counts.

    This function takes two dictionaries containing fragment formulas as keys and their counts 
    as values, merging them into a single dictionary with combined counts.

    Parameters:
    -----------
    counts1 : dict
        The first dictionary of fragment formulas and counts.

    counts2 : dict
        The second dictionary of fragment formulas and counts.

    Returns:
    --------
    dict
        A dictionary containing the combined fragment formulas and their counts.

    Notes:
    ------
    - If a formula appears in both dictionaries, its counts are summed.
    - If a formula appears only in one dictionary, it is included in the output as is.
    """
    combined_counts = counts1.copy()
    for formula, count in counts2.items():
        if formula in combined_counts:
            combined_counts[formula] += count
        else:
            combined_counts[formula] = count
    return combined_counts

def write_fragment_counts(file_path, fragment_formulas):
    """
    Writes fragment formulas and their counts to a specified output file.

    This function takes a dictionary of fragment formulas and their counts, sorts them in 
    descending order based on the counts, and writes the sorted results to the specified file.

    Parameters:
    -----------
    file_path : str
        Path to the output file where the fragment formulas and counts will be written.

    fragment_formulas : dict
        A dictionary where keys are fragment formulas (str) and values are their counts (int).

    Notes:
    ------
    - The output file will contain a title line followed by the formulas and counts.
    - The formulas are sorted from most common to least common based on their counts.
    """
    # Sort fragment formulas based on counts in descending order
    sorted_formulas = sorted(fragment_formulas.items(), key=lambda x: x[1], reverse=True)

    with open(file_path, 'w') as file:
        file.write("Fragment Formulas (Most to Least Common):\n")
        for formula, count in sorted_formulas:
            file.write(f"{formula}: {count}\n")

import json
from collections import defaultdict

def write_graphs_json():

    """
    Generates a JSON configuration for bond graphs and saves it to a file.

    This function reads bond data and molecular geometry from specified files to create a 
    configuration that describes bond relationships in a molecular structure. It generates 
    output files for each bond type and saves a JSON configuration that specifies 
    the relationships and properties of the bonds.

    Process:
    ---------
    - Reads bond data from 'bondarr.txt' to categorize bonds by type.
    - Extracts molecular geometry from a file named 'Geometry'.
    - Combines bond data specific to carbon atoms and aggregates them into separate files.
    - Configures and organizes bond types for graph generation, assigning colors for visualization.
    - Saves the entire configuration as a JSON file for use in graph generation.

    Raises:
    -------
    ValueError: If the 'Geometry' file does not contain a 'momentum' line, indicating 
    that the file does not have the expected format.

    Notes:
    ------
    - Assumes the existence of specific directory structures and file formats.
    - The output files for bond types are written to a results directory specified in the code.
    """
    bond_counts = defaultdict(list)

    # Reading bondarr.txt from ../results/bonds/
    with open('../results/bondarr.txt', 'r') as file:
        for line in file:
            bond, bond_type = line.strip().split(':')
            bond_counts[bond_type].append(bond)

    file_path = "Geometry"
    with open(file_path, 'r') as file:
        lines = file.readlines()

    natoms = 0
    for i, line in enumerate(lines):
        if line.strip().lower() == 'momentum':
            natoms = i
            break

    if natoms == 0:
        raise ValueError("The file does not contain a 'momentum' line.")

    geometry_lines = lines[0:natoms]
    symbols = [line.split()[0] for line in geometry_lines]

    carbon_files = defaultdict(list)
    carbon_counts = defaultdict(int)

    for i, symbol in enumerate(symbols):
        if symbol == "C":
            carbon_index = i + 1
            relevant_bonds = [
                f"../results/bonds/{bond}.out" for bond_type, bonds in bond_counts.items()
                for bond in bonds
                if str(carbon_index) in bond.split('-')
            ]
            combined_file = f"../results/bonds/carbon{carbon_index}.out"
            with open(combined_file, 'w') as outfile:
                for bond_file in relevant_bonds:
                    try:
                        with open(bond_file, 'r') as infile:
                            outfile.write(infile.read())
                    except FileNotFoundError:
                        print(f"Warning: {bond_file} not found.")
            carbon_files[carbon_index] = combined_file
            carbon_counts[carbon_index] = len(relevant_bonds)

    config = {
        "All": {
            "output_file": "../results/graphs/All.png",
            "files": []
        },
        "Carbons": {
            "output_file": "../results/graphs/Carbons.png",
            "files": []
        }
    }

    color_map = {
        "C-O": "red",
        "O-C": "red",
        "C-H": "blue",
        "C-F": "green",
        "F-C": "green",
        "C-C": "black",
        "C=C": "orange"
    }

    carbon_colors = ["red", "blue", "green", "black", "purple", "orange", "brown", "cyan", "magenta"]

    # Generate the "All" graph
    for bond_type, bonds in bond_counts.items():
        bond_filename = f"../results/bonds/{bond_type}.out"
        with open(bond_filename, 'w') as outfile:
            for bond in bonds:
                bond_file = f"../results/bonds/{bond}.out"
                try:
                    with open(bond_file, 'r') as infile:
                        outfile.write(infile.read())
                except FileNotFoundError:
                    print(f"Warning: {bond_file} not found.")
        config["All"]["files"].append({
            "filename": f"../results//bonds/{bond_type}.out",
            "label": bond_type,
            "no_bonds": len(bonds),
            "color": color_map.get(bond_type, "gray")
        })

    # Generate the "Carbons" graph
    for index, (carbon_index, combined_file) in enumerate(carbon_files.items()):
        config["Carbons"]["files"].append({
            "filename": combined_file,
            "label": f"Carbon {carbon_index}",
            "no_bonds": carbon_counts[carbon_index],
            "color": carbon_colors[index % len(carbon_colors)]
        })

    # Generate the bond-type-specific graphs
    for bond_type, bonds in bond_counts.items():
        config[bond_type] = {
            "output_file": f"../results/graphs/{bond_type}.png",
            "files": []
        }

        # Group bonds by the carbon atom in the bond
        carbon_env_files = defaultdict(list)
        for bond in bonds:
            atoms = bond.split('-')
            carbon_atoms = [atom for atom in atoms if symbols[int(atom) - 1] == 'C']
            for carbon in carbon_atoms:
                carbon_env_files[carbon].append(bond)

        # Combine files and add to the config
        for carbon, bond_list in carbon_env_files.items():
            combined_file = f"../results/bonds/{bond_type}_{carbon}.out"
            with open(combined_file, 'w') as outfile:
                for bond in bond_list:
                    bond_file = f"../results/bonds/{bond}.out"
                    try:
                        with open(bond_file, 'r') as infile:
                            outfile.write(infile.read())
                    except FileNotFoundError:
                        print(f"Warning: {bond_file} not found.")
            config[bond_type]["files"].append({
                "filename": f"../results/bonds/{bond_type}_{carbon}.out",
                "label": f"{bond_type} ({carbon})",
                "no_bonds": len(bond_list),
                "color": carbon_colors[int(carbon) % len(carbon_colors)]
            })

    # Save the JSON configuration to a file
    with open('../results/graphs_config.json', 'w') as json_file:
        json.dump(config, json_file, indent=4)

def track_completed_trajectory():
    """
    Tracks the number of completed trajectories by incrementing a count in a text file.

    This function updates a tracking file that maintains a count of how many trajectories 
    have been completed. It reads the current count from a file, increments it, and writes 
    the updated count back to the file.

    Process:
    ---------
    - If the tracking file does not exist, it initializes the count to zero.
    - Reads the current count from the file.
    - Increments the count and updates the file with the new count.

    Notes:
    ------
    - The tracking file is named 'completed_trajectories.txt' and is expected to be located 
      in the '../results/' directory.
    - The function prints the updated count of completed trajectories for user feedback.
    """
    # Path to the tracking file
    tracking_file = '../results/completed_trajectories.txt'
    
    # Initialize count if file doesn't exist
    if not os.path.exists(tracking_file):
        with open(tracking_file, 'w') as f:
            f.write("0")
    
    # Read the current count
    with open(tracking_file, 'r') as f:
        count = int(f.read().strip())
    
    # Increment the count
    count += 1
    
    # Write the updated count back to the file
    with open(tracking_file, 'w') as f:
        f.write(str(count))
    
    print(f"Trajectory completion recorded. Total completed: {count}")
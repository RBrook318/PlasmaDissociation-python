# Result file for the pure python plasma code. Should use a provided bondingarray file to know which bonds to investigate. Then it can look through xyz.all and 
# then make a dissociation.out as well as the C-C.out order files both in the output folder and in the total folder of the a whole set of trajectories. 
import os 

def process_results():
    # 1.read in the bonding array
    bondarr = read_bondarr()
    # # 2. detect dissociation
    detect_dissociation(bondarr)
    # 3 add to compiled results.
    compile_results()


def read_bondarr():
    bondarr = {}
    with open('../../results/bondarr.txt', 'r') as file:
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

    # Read the input data from a file
    with open('../output/xyz.all', 'r') as f:
        lines = f.readlines()

    # Initialize variables to store the current timestep, atom data, and dissociation flag
    timestep = None
    atoms = {}
    dissociated_bonds = set()  # Keep track of dissociated bond pairs

    # Open the output file for writing
    with open('../output/dissociation.out', 'w') as output:
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
    with open('../output/dissociation.out', "r") as f:
        lines = f.readlines()
        print(lines)

    for line in lines:
        if line.startswith("Dissociation detected"):
            parts = line.split(", ")
            timestep = int(float(parts[0].split(" ")[-1]))
            bond_info = parts[1].split(":")
            bond_type = bond_info[2].strip() # Remove leading and trailing whitespace
            bond_number = bond_info[1].strip()
            
            # Write to bond-type-specific output file
            bond_type_file = os.path.join("../../results", f"{bond_type}.out")
            with open(bond_type_file, "a") as f_out:
                f_out.write(f"{timestep}\n")
                

            
            # Write to old output file format and order by timestep
            bond_number_file = os.path.join("../../results", f"{bond_number}.out")
            with open(bond_number_file, "a") as f_out:
                f_out.write(f"{timestep}\n")

            
        



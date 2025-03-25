import re
from collections import defaultdict

def read_bond_types(file_path):
    bond_types = {}
    with open(file_path, 'r') as file:
        for line in file:
            line = line.strip()
            if line:
                bond_id, bond_type = line.split(':')
                bond_types[bond_id.strip()] = bond_type.strip()
    return bond_types

def parse_trajectory_data(file_path, bond_types):
    bond_breaks = defaultdict(list)
    current_trajectory = []
    trajectory_count = 0

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('--- run-'):
                if current_trajectory:
                    process_trajectory(current_trajectory, bond_types, bond_breaks, trajectory_count)
                    current_trajectory = []
                    trajectory_count += 1
                current_trajectory.append(line.strip())
            else:
                current_trajectory.append(line.strip())

        if current_trajectory:
            process_trajectory(current_trajectory, bond_types, bond_breaks, trajectory_count)

    return bond_breaks

def process_trajectory(trajectory, bond_types, bond_breaks, trajectory_count):
    if len(trajectory) == 0 or not trajectory[0].startswith('--- run-'):
        return

    trajectory_name = f"Trajectory {trajectory_count + 1}"

    for line in trajectory[1:]:
        match = re.search(r'Dissociation detected at timestep (\d+\.\d+), Broken bond: ([^:]+):([^\s]+)', line)
        if match:
            timestep = float(match.group(1))
            bond = match.group(2).strip()
            bond_type = bond_types.get(bond, 'Unknown')
            bond_breaks[trajectory_name].append((timestep, bond, bond_type))

def print_bond_breaks(bond_breaks):
    for trajectory, breaks in bond_breaks.items():
        print(f"Trajectory: {trajectory}")
        for timestep, bond, bond_type in breaks:
            print(f"    Timestep: {timestep}, Bond: {bond}, Bond Type: {bond_type}")
        print("\n")

def first_bond(bond_breaks): 
    first_type = defaultdict(list)
    first_number = defaultdict(list)

    for trajectory, breaks in bond_breaks.items():
        if breaks:
            timestep, bond, bond_type = breaks[0]
            first_type[bond_type].append(timestep)
            first_number[bond].append(timestep)

    total_first_bonds = sum(len(timesteps) for timesteps in first_type.values())

    with open('../../results/specifics/firstbond.out', 'w') as file:
        file.write(f"Type of bond most often broken first: \n\n")
        for bond_type, timesteps in first_type.items():
            percentage = (len(timesteps) / total_first_bonds) * 100
            avg_timestep = sum(timesteps) / len(timesteps)
            file.write(f"Bond type : {bond_type}, number of first breakages: {len(timesteps)} ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
        
        file.write(f"\nSpecific bond most often broken first: \n\n")
        for bond_number, timesteps in first_number.items():
            percentage = (len(timesteps) / total_first_bonds) * 100
            avg_timestep = sum(timesteps) / len(timesteps)
            file.write(f"Bond : {bond_number}, number of first breakages: {len(timesteps)} ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")

def second_bond(bond_breaks): 
    second_type = defaultdict(list)
    second_number = defaultdict(list)

    for trajectory, breaks in bond_breaks.items():
        if len(breaks) > 1:
            timestep, bond, bond_type = breaks[1]
            second_type[bond_type].append(timestep)
            second_number[bond].append(timestep)

    total_second_bonds = sum(len(timesteps) for timesteps in second_type.values())

    with open('../../results/specifics/secondbond.out', 'w') as file:
        file.write(f"Type of bond most often broken second: \n\n")
        for bond_type, timesteps in second_type.items():
            percentage = (len(timesteps) / total_second_bonds) * 100
            avg_timestep = sum(timesteps) / len(timesteps)
            file.write(f"Bond type : {bond_type}, number of second breakages: {len(timesteps)} ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
        
        file.write(f"\nSpecific bond most often broken second: \n\n")
        for bond_number, timesteps in second_number.items():
            percentage = (len(timesteps) / total_second_bonds) * 100
            avg_timestep = sum(timesteps) / len(timesteps)
            file.write(f"Bond : {bond_number}, number of second breakages: {len(timesteps)} ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")

def subsequent_bonds(bond_breaks, bond_types):
    subsequent_type_counts = defaultdict(int)
    subsequent_number_counts = defaultdict(int)
    bond_timesteps = defaultdict(list)

    for trajectory, breaks in bond_breaks.items():
        seen_bonds = set()
        
        for timestep, bond, bond_type in breaks:
            seen_bonds.add(bond)
            
            for seen_bond in seen_bonds:
                if seen_bond != bond:
                    seen_bond_type = bond_types.get(seen_bond, 'Unknown')
                    bond_type = bond_types.get(bond, 'Unknown')
                    subsequent_type_counts[(seen_bond_type, bond_type)] += 1
                    subsequent_number_counts[(seen_bond, bond)] += 1
                    bond_timesteps[(seen_bond, bond)].append(timestep)

    # Writing results to a file
    with open('../../results/specifics/subsequent_bonds_any.out', 'w') as file:
        file.write(f"Subsequent bond type counts (any subsequent bond):\n")
        for (bond1_type, bond2_type), count in subsequent_type_counts.items():
            total_bond1 = sum(v for k, v in subsequent_type_counts.items() if k[0] == bond1_type)
            percentage = (count / total_bond1) * 100
            if len(bond_timesteps[(bond1_type, bond2_type)]) > 0:
                avg_timestep = sum(bond_timesteps[(bond1_type, bond2_type)]) / len(bond_timesteps[(bond1_type, bond2_type)])
                file.write(f"Bond type {bond1_type} -> {bond2_type}: {count} occurrences ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
            else:
                file.write(f"Bond type {bond1_type} -> {bond2_type}: {count} occurrences ({percentage:.2f}%, avg. timestep = N/A)\n")

        file.write(f"\nSubsequent bond number counts (any subsequent bond):\n")
        for (bond1, bond2), count in subsequent_number_counts.items():
            total_bond1 = sum(v for k, v in subsequent_number_counts.items() if k[0] == bond1)
            percentage = (count / total_bond1) * 100
            if len(bond_timesteps[(bond1, bond2)]) > 0:
                avg_timestep = sum(bond_timesteps[(bond1, bond2)]) / len(bond_timesteps[(bond1, bond2)])
                file.write(f"Bond {bond1} -> {bond2}: {count} occurrences ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
            else:
                file.write(f"Bond {bond1} -> {bond2}: {count} occurrences ({percentage:.2f}%, avg. timestep = N/A)\n")

def next_bond(bond_breaks, bond_types):
    subsequent_type_counts = defaultdict(int)
    subsequent_number_counts = defaultdict(int)
    bond_timesteps = defaultdict(list)
    first_bond_count = defaultdict(int)

    for trajectory, breaks in bond_breaks.items():
        if len(breaks) > 1:
            for i in range(len(breaks) - 1):
                bond1 = breaks[i][1]
                bond2 = breaks[i + 1][1]
                timestep2 = breaks[i + 1][0]
                
                bond1_type = bond_types.get(bond1, 'Unknown')
                bond2_type = bond_types.get(bond2, 'Unknown')

                subsequent_type_counts[(bond1_type, bond2_type)] += 1
                subsequent_number_counts[(bond1, bond2)] += 1
                bond_timesteps[(bond1, bond2)].append(timestep2)
                first_bond_count[bond1] += 1

    # Writing results to a file
    with open('../../results/specifics/next_bonds.out', 'w') as file:
        file.write("Subsequent bond type counts:\n")
        for (bond1_type, bond2_type), count in subsequent_type_counts.items():
            total_bond1 = sum(v for k, v in subsequent_type_counts.items() if k[0] == bond1_type)
            percentage = (count / total_bond1) * 100
            if len(bond_timesteps[(bond1_type, bond2_type)]) > 0:
                avg_timestep = sum(bond_timesteps[(bond1_type, bond2_type)]) / len(bond_timesteps[(bond1_type, bond2_type)])
                file.write(f"Bond type {bond1_type} -> {bond2_type}: {count} occurrences ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
            else:
                file.write(f"Bond type {bond1_type} -> {bond2_type}: {count} occurrences ({percentage:.2f}%, avg. timestep = N/A)\n")

        file.write("\nSubsequent bond number counts:\n")
        sorted_keys = sorted(subsequent_number_counts.keys())
        for bond1, bond2 in sorted_keys:
            count = subsequent_number_counts[(bond1, bond2)]
            total_bond1 = sum(v for k, v in subsequent_number_counts.items() if k[0] == bond1)
            percentage = (count / first_bond_count[bond1]) * 100
            if len(bond_timesteps[(bond1, bond2)]) > 0:
                avg_timestep = sum(bond_timesteps[(bond1, bond2)]) / len(bond_timesteps[(bond1, bond2)])
                file.write(f"Bond {bond1} -> {bond2}: {count} occurrences ({percentage:.2f}%, avg. time = {avg_timestep*0.02418884254:.2f})\n")
            else:
                file.write(f"Bond {bond1} -> {bond2}: {count} occurrences ({percentage:.2f}%, avg. timestep = N/A)\n")

# Example usage
def specifics(bond_types_file):
    trajectory_data_file = '../../results/collated_diss.txt'
    bond_types = read_bond_types(bond_types_file)
    bond_breaks = parse_trajectory_data(trajectory_data_file, bond_types)
    # print_bond_breaks(bond_breaks)
    first_bond(bond_breaks)
    second_bond(bond_breaks)
    subsequent_bonds(bond_breaks, bond_types)
    next_bond(bond_breaks, bond_types)

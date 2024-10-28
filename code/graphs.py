import matplotlib.pyplot as plt
import json

def plot_data(filename, label,  color, total_trajectories):
    # Read data from file
    with open(filename, 'r') as file:
        lines = file.readlines()

    # Extract timesteps and sort them
    timesteps = [int(float(line.strip())) * 0.02418884254 for line in lines]
    timesteps.sort()

    # Calculate cumulative number of bonds broken
    cumulative_bonds = [0] * len(timesteps)
    for i, timestep in enumerate(timesteps):
        cumulative_bonds[i] = cumulative_bonds[i - 1] + 1 if i > 0 else 1

    # Calculate average number of bonds broken
    avg_bonds = [count / (total_trajectories) for count in cumulative_bonds]

    # Plot average number of bonds broken
    plt.plot(timesteps, avg_bonds, label=label, color=color)

def generate_graphs_from_json(json_file):
    # Load the JSON file
    with open(json_file, 'r') as file:
        graphs = json.load(file)
    
    # Read the total number of bonds from allbonds.out
    with open('../results/completed_trajectories.txt', 'r') as all_file:
        lines = all_file.readlines() 
        total_trajectories = lines  # Total number of lines in the file

    # Iterate through each graph configuration in the JSON
    for graph_name, graph_data in graphs.items():
        plt.figure()

        # Iterate through each file info in the current graph data
        for file_info in graph_data['files']:
            # Adjust the path to look for the files in ../results/bonds/
            bond_file_path = file_info['filename']
            plot_data(bond_file_path, file_info['label'], 
                      file_info['color'], total_trajectories)

        # Set the graph labels and title
        plt.xlabel('Femtoseconds (fs)')
        plt.ylabel('Percentage of bonds broken (%)')
        plt.title(f"Bonds broken for {graph_name}")
        plt.legend()

        # Save the plot to ../results/graphs/
        output_path = graph_data['output_file']
        plt.savefig(output_path)
        plt.close()

def create_graphs():
    json_file = '../results/graphs_config.json'  # Path to your JSON file
    generate_graphs_from_json(json_file)



import matplotlib.pyplot as plt
import json

def plot_data(filename, label, no_bonds, color, total_trajectories):
    """
    Plots the average number of bonds broken per trajectory from data in a specified file.

    This function reads bond data from the given file, calculates the cumulative number of 
    bonds broken, and then plots the average number of bonds broken as a function of time.

    Parameters:
    -----------
    filename : str
        The path to the file containing bond data, with each line representing a timestep.
    label : str
        The label for the plot line, used in the legend.
    color : str
        The color for the plot line.
    total_trajectories : int
        The total number of trajectories, used to calculate the average number of bonds broken.

    Notes:
    ------
    - The function expects the file to contain lines that can be converted to float values.
    - The time unit is converted to femtoseconds using a specific conversion factor.
    """
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
    percentage_bonds = [100*count/(total_trajectories*no_bonds) for count in cumulative_bonds]

    # Plot average number of bonds broken
    plt.plot(timesteps, percentage_bonds, label=label, color=color)

def generate_graphs_from_json(json_file):
    """
    Generates bond graphs based on configurations specified in a JSON file.

    This function reads a JSON configuration file to retrieve graph details, including 
    file paths and settings for each bond type. It then calls `plot_data` for each file 
    to generate and save the corresponding plots.

    Parameters:
    -----------
    json_file : str
        The path to the JSON file containing graph configurations.

    Notes:
    ------
    - The function reads from '../results/completed_trajectories.txt' to obtain the total 
      number of trajectories.
    - It expects the JSON structure to contain specific fields for each graph.
    - The generated plots are saved in the specified output paths defined in the JSON file.
    """
    # Load the JSON file
    with open(json_file, 'r') as file:
        graphs = json.load(file)
    
    # Read the total number of bonds from allbonds.out
    with open('../results/completed_trajectories.txt', 'r') as all_file:
        total_trajectories = int(all_file.read().strip()) 

    # Iterate through each graph configuration in the JSON
    for graph_name, graph_data in graphs.items():
        plt.figure()

        # Iterate through each file info in the current graph data
        for file_info in graph_data['files']:
            # Adjust the path to look for the files in ../results/bonds/
            bond_file_path = file_info['filename']
            plot_data(bond_file_path, file_info['label'], file_info['no_bonds'],
                      file_info['color'], total_trajectories)

        # Set the graph labels and title
        plt.xlabel('Femtoseconds (fs)')
        plt.ylabel('Average number of bonds broken per trajectory')
        plt.title(f"Bonds broken for {graph_name}")
        plt.legend()

        # Save the plot to ../results/graphs/
        output_path = graph_data['output_file']
        plt.savefig(output_path)
        plt.close()

def create_graphs():
    """
    Creates and saves bond graphs based on the configurations in the JSON file.

    This function serves as the entry point to generate graphs. It specifies the path to 
    the JSON configuration file and calls `generate_graphs_from_json` to produce the 
    graphs.

    Notes:
    ------
    - The JSON file is expected to be located at '../results/graphs_config.json'.
    - This function does not return any values but generates and saves plots.
    """
    json_file = '../results/graphs_config.json'  # Path to your JSON file
    generate_graphs_from_json(json_file)

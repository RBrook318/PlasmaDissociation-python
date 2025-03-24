import json

# Initialize global variables
cores = None
num_states = None
multiplicity = None
start_state = None
method = None
spin_flip = None
basis = None
gpu = None
branches = None
timestep = None
tot_timesteps = None
checks = None
remove_atoms = None

# Function to load variables from JSON file
def load_global_variables():
    global cores,num_states, multiplicity, start_state
    global method, spin_flip, basis, gpu
    global branches, timestep, tot_timesteps, checks, remove_atoms

    # Open and read the JSON file
    with open("../../inputs.json", "r") as f:
        inputs = json.load(f)
    
    # Set global variables from JSON file
    cores = inputs["HPC"]["cores"]
    
    num_states = inputs["Molecule_data"]["States"]
    multiplicity = inputs["Molecule_data"]["Multiplicity"]
    start_state = inputs["Molecule_data"]["Start_State"]

    
    method = inputs["Elec_structure"]["method"]
    spin_flip = inputs["Elec_structure"]["Spin_flip"]
    basis = inputs["Elec_structure"]["Basis"]
    gpu = inputs["Elec_structure"]["GPU"]
    
    branches = inputs["Propagation"]["Branches"]
    timestep = inputs["Propagation"]["Timestep"]
    tot_timesteps = inputs["Propagation"]["Tot_timesteps"]
    checks = inputs["Propagation"]["Checks"]
    remove_atoms = inputs["Propagation"]["Remove_atoms"]
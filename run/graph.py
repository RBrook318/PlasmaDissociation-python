
import os
import matplotlib.pyplot as plt
import matplotlib as mpl
import json 

with open('inputs.json') as f:
        inputs=json.load(f)
title=inputs["run"]["Molecule"]
reps=float(inputs["setup"]["repeats"])
colors =mpl.colormaps['prism']
def plot_data(filename_in, filename_out, label, no_bonds):
    # Read data from file
    with open(filename_in, 'r') as file:
        lines = file.readlines()

    # Extract timesteps and sort them
    timesteps = [int(line.strip()) for line in lines]
    timesteps.sort()

    # Calculate cumulative number of bonds broken
    cumulative_bonds = [0] * len(timesteps)
    for i, timestep in enumerate(timesteps):
        cumulative_bonds[i] = cumulative_bonds[i-1] + 1 if i > 0 else 1

    # Calculate average number of bonds broken
    avg_bonds = [count / no_bonds for count in cumulative_bonds]
    scaled_avg_bonds = [value * 100 for value in avg_bonds]
    # Plot average number of bonds broken
    plt.plot(timesteps, scaled_avg_bonds, label=label)
    plt.xlabel('Timesteps')
    plt.ylabel('Bonds broken')
    plt.title('Number of ' + label +' bonds broken in '+title)
    plt.ylim(0, 100)  # Set the y-axis limits
    plt.xlim(0)  # Set the x-axis lower limit to 0
    plt.legend()
    # Display the plot
    plt.savefig(filename_out)
    plt.close()

def plot_data_cum(folder_path, file_contents, filename_in, filename_out, label, no_bonds):
    # Read data from file
    with open(filename_in, 'r') as file:
        lines = file.readlines()

    # Extract timesteps and sort them
    timesteps = [int(line.strip()) for line in lines]
    timesteps.sort()
    # Calculate cumulative number of bonds broken
    cumulative_bonds = [0] * len(timesteps)
    for i, timestep in enumerate(timesteps):
        cumulative_bonds[i] = cumulative_bonds[i-1] + 1 if i > 0 else 1
    # Calculate average number of bonds broken
    avg_bonds = [count / no_bonds for count in cumulative_bonds]
    scaled_avg_bonds_all = [value * 100 for value in avg_bonds]
    # Plot average number of bonds broken
    plt.plot(timesteps, scaled_avg_bonds_all, label=label+' average')
    leg_cnt=0
    style=['dashed','dotted','dashdot']
    style_choice=style[0]
    alpha=0.8
    for name in file_contents.split('\n'):
        if name.endswith(label): 
            leg_cnt=leg_cnt+1
            if(alpha<0.3):
                alpha=0.8
            if(leg_cnt>10):
                style_choice=style[1]
            if(leg_cnt>20):
                style_choice=style[2]
            filename=name[0:len(name)-len(label)-1]
            with open(folder_path+filename+'.out', 'r') as file:
                lines = file.readlines()
                # Extract timesteps and sort them
                timesteps = [int(line.strip()) for line in lines]
                timesteps.sort()
                # Calculate cumulative number of bonds broken
                cumulative_bonds = [0] * len(timesteps)
                for i, timestep in enumerate(timesteps):
                    cumulative_bonds[i] = cumulative_bonds[i-1] + 1 if i > 0 else 1
                # Calculate average number of bonds broken
                avg_bonds = [count/reps for count in cumulative_bonds]
                scaled_avg_bonds = [value * 100 for value in avg_bonds]
                # Plot average number of bonds broken
                plt.plot(timesteps, scaled_avg_bonds, label=filename, alpha=0.5, linestyle=style_choice,color=colors(leg_cnt*10))
                alpha=alpha-0.1                                  
    plt.xlabel('Timesteps')
    plt.ylabel('Bonds broken')
    plt.title('Number of ' + label +' bonds broken in '+title)
    if leg_cnt>10:
        plt.legend(ncol=2)
    else:
        plt.legend()
    plt.ylim(0, 100)  # Set the y-axis limits
    plt.xlim(0)
    # Display the plot
    plt.savefig(filename_out)
    plt.close()

def plot_data_all(folder_path,file_contents,filename_out):
    for filename in os.listdir(folder_path):
        if os.path.isfile(os.path.join(folder_path, filename)):
            if filename.endswith('.out'):
                filename_without_extension = filename[:-4]
                if filename[0].isalpha():
                    bond_cnt=0
                    for line in file_contents.split('\n'):
                        if line.endswith(filename_without_extension):
                            bond_cnt+=1
                    with open(folder_path+filename, 'r') as file:
                        lines = file.readlines()
                    bond_cnt=bond_cnt*reps
                    # Extract timesteps and sort them
                    timesteps = [int(line.strip()) for line in lines]
                    timesteps.sort()
                    # Calculate cumulative number of bonds broken
                    cumulative_bonds = [0] * len(timesteps)
                    for i, timestep in enumerate(timesteps):
                        cumulative_bonds[i] = cumulative_bonds[i-1] + 1 if i > 0 else 1

                    # Calculate average number of bonds broken
                    avg_bonds = [count / bond_cnt for count in cumulative_bonds]
                    scaled_avg_bonds = [value * 100 for value in avg_bonds]
                    # Plot average number of bonds broken
                    plt.plot(timesteps, scaled_avg_bonds, label=filename_without_extension)
    plt.xlabel('Timesteps')
    plt.ylabel('Bonds broken')
    plt.ylim(0, 100)  # Set the y-axis limits
    plt.xlim(0)
    plt.legend()
    plt.savefig(filename_out)
    plt.close()

with open('results/bondarr.txt', 'r') as file:
    file_contents = file.read()

folder_path = 'results/'  # Replace with the actual folder path
os.mkdir('results/graphs')
plot_folder='results/graphs/'
# Iterate over the files in the folder
for filename in os.listdir(folder_path):
    if os.path.isfile(os.path.join(folder_path, filename)):
        if filename.endswith('.out'):
            filename_without_extension = filename[:-4]
            if filename[0].isalpha():
                bond_cnt=0
                for line in file_contents.split('\n'):
                    if line.endswith(filename_without_extension):
                        bond_cnt+=1
                # Plot data from filenames starting with a letter
                plot_data(os.path.join(folder_path, filename),plot_folder+filename_without_extension+'.png', filename_without_extension, reps*bond_cnt)
                
            else:
                # Find line in file_contents that starts with file_without_extension
                for line in file_contents.split('\n'):
                    if line.startswith(filename_without_extension):
                        full_name=line
                plot_data(os.path.join(folder_path, filename),plot_folder+full_name+'.png', full_name, reps)
                # Plot data from filenames starting with a number


for filename in os.listdir(folder_path):
    if os.path.isfile(os.path.join(folder_path, filename)):
        if filename.endswith('.out'):
            filename_without_extension = filename[:-4]
            if filename[0].isalpha():
                bond_cnt=0
                for line in file_contents.split('\n'):
                    if line.endswith(filename_without_extension):
                        bond_cnt+=1
                plot_data_cum(folder_path,file_contents, os.path.join(folder_path, filename),plot_folder+filename_without_extension+'_cum.png', filename_without_extension, reps*bond_cnt)


plot_data_all(folder_path,file_contents,plot_folder+'all.png')

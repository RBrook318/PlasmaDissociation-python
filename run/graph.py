import matplotlib.pyplot as plt

def plot_data(filename, label, no_bonds):
    # Read data from file
    with open(filename, 'r') as file:
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

    # Plot average number of bonds broken
    plt.plot(timesteps, avg_bonds, label=label)


# Plot data from 'case1.out'
plot_data('C-H.out', 'C-H', 3)

# Plot data from 'Case1c500.out'
# plot_data('2-4.out', 'Lower H', 1)

# Plot data from 'Case1c1000.out'
# plot_data('C-H.out', 'C-H', 2)

# Add labels and title to the plot
plt.xlabel('Timesteps')
plt.ylabel('Bonds broken')
plt.title('Number of each hydrogen bond broken in C3H2F6')
plt.legend()

# Display the plot
plt.savefig('C3H2F6-Hbonds.png')

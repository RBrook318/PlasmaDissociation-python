# Define the mass of an electron in atomic units (assuming the mass unit is 1.0)
mass = [12, 19, 19, 19, 12, 19, 12, 19, 19]
for i in range(len(mass)):
    mass[i] = mass[i] * 1836

# Initialize variables to store total kinetic energy and the number of files processed
total_kinetic_energy = 0.0
num_files = 500
atoms = 9

# Loop through all the geometry files (from 1 to 500)
for x in range(1, num_files + 1):
    # Read the momentum values from the current file (assuming they are in the provided format)
    file_path = f'Geometry.{x}'
    momentum_lines = open(file_path).readlines()[-atoms:-1]
    
    n = 0

    # Initialize variables to store kinetic energy for the current file
    file_kinetic_energy = 0.0

    # Loop through each line of momentum data
    for line in momentum_lines:
        # Split the line into components and convert them to floating-point numbers
        components = [float(value) for value in line.split()]

        # Calculate the magnitude of the momentum vector (sqrt(p_x^2 + p_y^2 + p_z^2))
        magnitude = (components[0]**2 + components[1]**2 + components[2]**2)**0.5

        # Calculate the kinetic energy for this particle
        kinetic_energy = 0.5 * (magnitude**2 / mass[n])

        # Add the kinetic energy to the total for this file
        file_kinetic_energy += kinetic_energy
        n = n + 1

    # Calculate the average kinetic energy for this file and add it to the total 
    total_kinetic_energy += file_kinetic_energy

# Calculate the overall average kinetic energy
average_total_kinetic_energy = total_kinetic_energy / num_files

# Print the average total kinetic energy
print("Average Total Kinetic Energy: {:.6f} Hartrees".format(average_total_kinetic_energy))

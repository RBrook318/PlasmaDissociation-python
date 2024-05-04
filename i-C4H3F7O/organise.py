import numpy as np

# Load your data from the file
data = np.loadtxt('modes.txt', dtype=str)  # Load data as strings

# Extract only the numeric columns (ignoring the symbol column)
numeric_data = data[:, 1:].astype(float)

# Combine the columns
column1 = np.concatenate((numeric_data[:, 0], numeric_data[:, 3], numeric_data[:, 6]))
column2 = np.concatenate((numeric_data[:, 1], numeric_data[:, 4], numeric_data[:, 7]))
column3 = np.concatenate((numeric_data[:, 2], numeric_data[:, 5], numeric_data[:, 8]))

# Stack the columns horizontally to create a new array
reshaped_data = np.column_stack((column1, column2, column3))

# Save the reshaped data to a new file
np.savetxt('a.txt', reshaped_data, fmt='%f', delimiter='\t')

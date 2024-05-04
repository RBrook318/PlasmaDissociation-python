import numpy as np

# Define the coordinates of the atoms in the new geometry
coordinates = np.array([
    [-1.7363594157, 0.2935369014, 0.0148210222],
    [-2.4296590856, 0.1623666083, -1.1141104457],
    [-1.4900018788, 1.5823120339, 0.2059076529],
    [-2.5156754691, -0.1342967001, 1.0060215187],
    [-0.4771558046, -0.5052274620, -0.0487374809],
    [-0.6883189256, -1.8039897125, -0.2526449093],
    [0.7608778828, -0.0673949107, 0.0696511892],
    [1.1074519775, 1.1699518576, 0.2717736887],
    [1.7930692137, -0.8573805048, -0.0066829793]
])

# Define the bonding array (a list of tuples indicating which atoms are bonded)
bonding_array = [(0, 4), (0, 1), (0, 2), (0, 3), (4, 5), (4, 6), (6, 7), (6, 8)]

# Calculate and print bond distances for the new geometry
for bond in bonding_array:
    atom1 = coordinates[bond[0]]
    atom2 = coordinates[bond[1]]
    distance = np.linalg.norm(atom1 - atom2)
    print(f'Bond distance between atoms {bond[0] + 1} and {bond[1] + 1}: {distance:.4f} Ångstroms')

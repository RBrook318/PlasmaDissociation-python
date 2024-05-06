# Function to read coordinates from a file
def read_coordinates_from_file(file_path):
    coordinates = []
    with open(file_path, 'r') as file:
        for line in file:
            if line.strip():  # Check if the line is not empty
                atom, x, y, z = line.split()
                coordinates.append((atom, float(x), float(y), float(z)))
    return coordinates

# Convert coordinates from angstrom to Bohr

def convert_to_bohr(coordinates):
    angstrom_to_bohr = 1.88973
    return [(atom, x * angstrom_to_bohr, y * angstrom_to_bohr, z * angstrom_to_bohr) for atom, x, y, z in coordinates]

# Read coordinates from file
file_path = 'geom.txt'
coordinates = read_coordinates_from_file(file_path)

# Convert coordinates to Bohr
coordinates_bohr = convert_to_bohr(coordinates)

# Display converted coordinates
for atom, x, y, z in coordinates_bohr:
    print(f"{atom:<2} {x:>12.10f} {y:>12.10f} {z:>12.10f}")
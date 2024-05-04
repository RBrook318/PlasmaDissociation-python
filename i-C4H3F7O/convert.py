# Conversion factor: 1 Bohr = 0.529177 angstroms
angstrom_to_bohr = 0.529177

# Function to convert coordinates from angstroms to Bohr
def convert_to_bohr(coordinates):
    bohr_coordinates = []
    for line in coordinates:
        elements = line.split()
        element = elements[0]
        x = float(elements[1]) / angstrom_to_bohr
        y = float(elements[2]) / angstrom_to_bohr
        z = float(elements[3]) / angstrom_to_bohr
        bohr_coordinates.append(f"{element:2s}  {x:.10f}  {y:.10f}  {z:.10f}")

    return bohr_coordinates

# Example usage
angstrom_coordinates = [
    "C  1.5838263664  -2.1138974121  -0.2376754673",
    "O  0.3239491845  -1.4735602944  -0.4400182213",
    "C  0.0268048348  -0.3512092542   0.2388887853",
    "C  0.9935906731   0.8242549235  -0.0595082763",
    "F  2.1693399812   0.5931588457   0.5173682063",
    "F  1.1919166258   0.9413300713  -1.3647985474",
    "F  0.5310313639   1.9754781668   0.4009850865",
    "C -1.4462406499  -0.0237836128  -0.0958410659",
    "F -2.1965396269  -1.0932554202   0.1051499807",
    "F -1.9026642495   0.9523069651   0.6739372807",
    "F -1.5649270224   0.3490875468  -1.3632842102",
    "F  0.0844312324  -0.5349303597   1.5911596710",
    "H  1.4636997212  -3.1034586800  -0.6559993066",
    "H  2.3710330405  -1.5847888740  -0.7628751242",
    "H  1.8124916749  -2.1840402996   0.8191791378",
]

bohr_coordinates = convert_to_bohr(angstrom_coordinates)

# Print the result
for line in bohr_coordinates:
    print(line)

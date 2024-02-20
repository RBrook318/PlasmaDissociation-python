import numpy as np

def calculate_difference(file1_path, file2_path, output_path):
    # Load coordinates from each file
    coordinates1 = np.loadtxt(file1_path)
    coordinates2 = np.loadtxt(file2_path)

    # Check if the shape of the arrays match
    if coordinates1.shape != coordinates2.shape:
        print("Error: Arrays have different shapes.")
        return

    # Calculate differences
    differences = abs(coordinates2 - coordinates1) 

    # Save differences to a new file
    np.savetxt(output_path, differences, fmt='%.20f')


file1_path = "python.txt"
file2_path = "fortran.txt"
output_path = "diffrences.txt"

calculate_difference(file1_path, file2_path, output_path)
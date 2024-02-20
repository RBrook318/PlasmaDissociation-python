import os

def process_geometry_file(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    # Remove the first line
    lines.pop(0)

    # Add amplitudes at the end of the file
    amplitudes = '\n 1.0 \n 0.0\n'
    lines.extend(amplitudes)

    # Write the modified lines back to the file
    with open(file_path, 'w') as file:
        file.writelines(lines)


# Iterate through all files in the directory
for i in range(1, 501):
    filename = f'Geometry.{i}'
    process_geometry_file(filename)

import os

# Directory where your geometry.x files are located
directory = os.getcwd()

# Loop through files from 1 to 500
for x in range(1, 501):
    # Construct the file path
    file_path = os.path.join(directory, f'geometry.{x}')

    # Check if the file exists
    if os.path.exists(file_path):
        # Read the contents of the file
        with open(file_path, 'r') as file:
            lines = file.readlines()

        # Check if the file has at least one line
        if len(lines) > 0:
            # Remove the empty first line
            lines = lines[1:]

            # Write the modified content back to the file
            with open(file_path, 'w') as file:
                file.writelines(lines)

        print(f'Removed empty first line from {file_path}')
    else:
        print(f'File {file_path} does not exist.')

print('Processing complete.')

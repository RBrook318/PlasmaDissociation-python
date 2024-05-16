#########################################################################################
#
#   Python Run script for setting up run paramters for dissociation after electron impact
#   Written by O Bramley                                          08.05.24
# 
#     
#
#########################################################################################


import numpy as np
import json
from pyqchem import QchemInput
from pyqchem.tools import get_geometry_from_pubchem
from pyqchem import get_output_from_qchem
from pyqchem.parsers.parser_optimization import basic_optimization
from pyqchem.parsers.parser_frequencies import basic_frequencies


def convert_to_bohr(coordinates):
    angstrom_to_bohr = 1.88973
    return [(atom, x * angstrom_to_bohr, y * angstrom_to_bohr, z * angstrom_to_bohr) for atom, x, y, z in coordinates]

def organise_modes(modes):
    numeric_modes = np.zeros((len(modes), len(modes[0])-1))
    for i in range(len(modes)):
        for j in range(len(modes[0])-1):
            numeric_modes[i,j] = float(modes[i][j+1])

    # Combine the columns
    column1 = np.concatenate((numeric_modes[:, 0], numeric_modes[:, 3], numeric_modes[:, 6]))
    column2 = np.concatenate((numeric_modes[:, 1], numeric_modes[:, 4], numeric_modes[:, 7]))
    column3 = np.concatenate((numeric_modes[:, 2], numeric_modes[:, 5], numeric_modes[:, 8]))
    # Stack the columns horizontally to create a new array
    reshaped_data = np.column_stack((column1, column2, column3))
    return reshaped_data

def bondarr(output):
    z_matrix = []

    for line in file:
        if line.strip() == "Z-matrix Print:":
            break
    file.readline()
    file.readline()
    for i in range(atoms):
        z_matrix.append(file.readline().split())

    bonds = np.zeros((atoms, atoms))
    for i in range(1,atoms):
        bonds[i,int(z_matrix[i][1])-1] = 1
    with open("../results/bondarr.txt", "w") as file:
        for i in range(atoms):
            for j in range(atoms):
                if bonds[j,i] == 1:
                    file.write((f'{i+1}'+'-'+f'{j+1}'+':'+z_matrix[i][0]+'-'+z_matrix[j][0]+'\n')) 

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:, 0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
    # Initialize arrays for random
    Meff = np.zeros(nmod)
    rv = np.zeros((nmod, mom_num))
    for i in range(nmod):
        for j in range(n):
            Meff[i] = Meff[i]+np.sum(((Ax[j, i]**2) + (Ay[j, i]**2) + (Az[j, i]**2)) * m[j])
        rv[i, :] = rn[i, :] * np.sqrt(2 * T / Meff[i])
    # Calculate the velocity by applying it through the tranformation matrix of normal modes.
    Vx = np.dot(Ax, rv)
    Vy = np.dot(Ay, rv)
    Vz = np.dot(Az, rv)
    Px = np.zeros((n,mom_num))
    Py = np.zeros((n,mom_num))
    Pz = np.zeros((n,mom_num))
    for i in range(n):
        Px[i,:] = Vx[i,:]*m[i]
        Py[i,:] = Vy[i,:]*m[i]
        Pz[i,:] = Vz[i,:]*m[i]
    
    return Px, Py, Pz

if __name__ == "__main__":
    with open('../inputs.json') as f:
        inputs=json.load(f)
    # Load molecule coordinates from pubchem
    molecule = get_geometry_from_pubchem(inputs["run"]["Molecule"])
    # Optimise the molecule
    qc_inp = QchemInput(molecule,
                        jobtype='opt',
                        exchange='BHHLYP !50% HF +  50% Becke88 exchange',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm='diis')
    
    output = get_output_from_qchem(qc_inp,processors=2)
    open('optimisation.out').write(output)
    pasrser_output = basic_optimization(output)
    opt_geoms=(pasrser_output['optimized_molecule'])
    # Write bond breaking file to results folder
    bondarr(output)

    # Find Normal modes
    qc_inp = QchemInput(opt_geoms,
                        jobtype='FREQ',
                        exchange='BHHLYP',
                        basis='6-31+G*')
    output = get_output_from_qchem(qc_inp,processors=2)
    open('modes.out').write(output)
    pasrser_output = basic_frequencies(output)
    num_modes = len(pasrser_output['number_of_modes'])
    # Need to check this works correctly 
    modes=organise_modes(pasrser_output['modes'])
    # Convert Geometries to bohr
    opt_geoms=convert_to_bohr(opt_geoms)
    open(inputs["run"]["Molecule"]+'.xyz').write(opt_geoms)

    # Extract masses of atoms   
    masses=int(qc_inp.molecule.get_atomic_numbers())
    Px, Py, Pz = create_geom(inputs["run"]["Atoms"],num_modes,inputs["run"]["Temp"],modes,masses,inputs["setup"]["repeats"])
    # Extract atom symbols
    atoms=qc_inp.molecule.get_symbols()
    # Write momenta files to repetition folder
    for j in range(inputs["setup"]["repeats"]):
        with open('../rep'+str(j+1)+'/Geometry', 'w') as file:
            for atom in range(inputs["run"]["Atoms"]):
               file.write(f'{atoms[atom]}  {opt_geoms[atom][0]}  {opt_geoms[atom][1]}  {opt_geoms[atom][2]}\n') 
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(inputs["run"]["Atoms"]):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')

# The geometry can by specified in bohr by setting the $rem variable INPUT_BOHR equal to TRUE.
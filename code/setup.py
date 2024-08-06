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
    values=coordinates.get_coordinates()
    atom=coordinates.get_symbols()
    output = ""
    for i in range(len(atom)):
        output+=(atom[i]+" "+str(float(values[i][0])*angstrom_to_bohr)+" "+str(float(values[i][1])*angstrom_to_bohr)+" "+str(float(values[i][2])*angstrom_to_bohr)+" \n")
    return output
   

def organise_modes(modes,atoms):
    numeric_modes = np.zeros((len(modes)*atoms, 3))
    cnt=0
    for i in range(len(modes)):
        for j in range(atoms):
            numeric_modes[cnt,0]= float(modes[i]['displacement'][j][0])
            numeric_modes[cnt,1]= float(modes[i]['displacement'][j][1])
            numeric_modes[cnt,2]= float(modes[i]['displacement'][j][2])
            cnt+=1
   
    return numeric_modes

def bondarr(molecule):
    bonds = molecule.get_connectivity()
    atoms = molecule.get_symbols()
    unique_bonds = set()
    
    with open("../results/bondarr.txt", "w") as file:
        for bond in bonds:
            sorted_bond = tuple(sorted(bond))
            if sorted_bond not in unique_bonds:
                unique_bonds.add(sorted_bond)
                file.write(f"{bond[0]}-{bond[1]}:{atoms[bond[0]-1]}-{atoms[bond[1]-1]}\n")
    
    return

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:,0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
    m=m*1822.8885300626
    T=T*0.0000031668
    print(T)
    # Initialize arrays for random
    Meff = np.zeros(nmod)
    rv = np.zeros((nmod, mom_num))
    for i in range(nmod):
        for j in range(n):
            print(m[j])
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
    # Write bond breaking file to results folder
    if inputs["run"]["method"] == "QChem":
        qc.initial_conditions
    # Optimise the molecule
    qc_inp = QchemInput(molecule,
                        jobtype='opt',
                        exchange='BHHLYP !50% HF +  50% Becke88 exchange',
                        basis='6-31+G*',
                        unrestricted=True,
                        max_scf_cycles=500,
                        sym_ignore=True,
                        scf_algorithm='diis')
    
    output = get_output_from_qchem(qc_inp,processors=inputs['setup']['cores'])
    with open('optimisation.out', 'w') as f:
        f.write(output)
    pasrser_output = basic_optimization(output)
    opt_geoms=(pasrser_output['optimized_molecule'])
    bondarr(opt_geoms)
    # Find Normal modes
    qc_inp = QchemInput(opt_geoms,
                        jobtype='FREQ',
                        exchange='BHHLYP',
                        basis='6-31+G*')
    output = get_output_from_qchem(qc_inp,inputs['setup']['cores'])
    with open('modes.out', 'w') as f:
        f.write(output)
    pasrser_output = basic_frequencies(output)
    num_modes = len(pasrser_output['modes'])
    # Need to check this works correctly 
    natoms = int((num_modes+6)/3)
    modes=organise_modes(pasrser_output['modes'],natoms)
    # Convert Geometries to bohr
    opt_geoms=convert_to_bohr(opt_geoms)
    with open(inputs["run"]["Molecule"]+'.xyz','w') as f:
        f.write(opt_geoms)

    # Extract masses of atoms   
    masses=(qc_inp.molecule.get_atomic_masses())
    Px, Py, Pz = create_geom(natoms,num_modes,inputs["run"]["Temp"],modes,masses,inputs["setup"]["repeats"])
    # Extract atom symbols
    atoms=qc_inp.molecule.get_symbols()
    # Write momenta files to repetition folder
    for j in range(inputs["setup"]["repeats"]):
        with open('../rep-'+str(j+1)+'/Geometry', 'w') as file:
            file.write(opt_geoms)
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(natoms):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')

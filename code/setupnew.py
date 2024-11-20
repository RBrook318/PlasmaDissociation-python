import numpy as np
import json
import pubchempy as pcp
import init
# from chemformula import ChemFormula

def create_geom(n,nmod,T,modes,m,mom_num):
    Ax = modes[:,0]
    Ay = modes[:,1]
    Az = modes[:,2]
    Ax = Ax.reshape(n, nmod, order = 'F')
    Ay = Ay.reshape(n, nmod, order = 'F')
    Az = Az.reshape(n, nmod, order = 'F')
    rn = np.random.randn(nmod, mom_num)  # Use np.random.randn for standard normal distribution
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

def get_geometry(molecule_name):
    search = pcp.get_compounds(molecule_name, 'name', record_type='3d')
    print(search)
    molecule = search[0]
    bonds = molecule.bonds
    atoms = molecule.atoms
    
    # Create a mapping from atom IDs to symbols
    atom_id_to_symbol = {atom.aid: atom.element for atom in atoms}
    symbols = [atom.element for atom in atoms]
    coordinates = [[atom.x*1.88973, atom.y*1.88973, atom.z*1.88973] for atom in atoms]
    with open("../results/bondarr.txt", "w") as file:
        for bond in bonds:
            atom1 = bond.aid1
            atom2 = bond.aid2
            order = bond.order
            
            # Determine the bond type based on the order
            if order == 1:
                bond_type = "-"
            elif order == 2:
                bond_type = "="
            elif order == 3:
                bond_type = "#"
            else:
                bond_type = "-"  # Default to single bond if unknown order
            
            # Get the symbols for the atoms
            atom1_symbol = atom_id_to_symbol.get(atom1, "Unknown")
            atom2_symbol = atom_id_to_symbol.get(atom2, "Unknown")
            
            # Write to file in the format 1-3:C-H
            file.write(f"{atom1}-{atom2}:{atom1_symbol}{bond_type}{atom2_symbol}\n")

    return symbols,coordinates

def organise_modes(modes,atoms):
    numeric_modes = np.zeros((len(modes)*atoms, 3))
    print(modes[0])
    cnt=0
    for i in range(len(modes)):
        for j in range(atoms):
            numeric_modes[cnt,0]= float(modes[i]['displacement'][j][0])
            numeric_modes[cnt,1]= float(modes[i]['displacement'][j][1])
            numeric_modes[cnt,2]= float(modes[i]['displacement'][j][2])
            cnt+=1
   
    return numeric_modes

if __name__ == "__main__":
    with open('../inputs.json') as f:
        inputs=json.load(f)
    symbols,coordinates = get_geometry(inputs["run"]["Molecule"])
    if inputs["run"]["method"] == "QChem":
        import qchem as qc
        opt_coords, modes = qc.initial_conditions(symbols,coordinates,inputs['setup']['cores'])
    elif inputs["run"]["method"] == "PySCF":
        import elec_pyscf as pyscf
        pyscf.initial_conditions(symbols,coordinates)
    
   
    modes = modes['modes']
    num_modes = len(modes)
    natoms = int((num_modes+6)/3)
    modes=organise_modes(modes,natoms)
    masses = init.setup_masses(symbols)
    Px, Py, Pz = create_geom(natoms,num_modes,inputs["run"]["Temp"],modes,masses,inputs["setup"]["repeats"])
    for j in range(inputs["setup"]["repeats"]):
        with open('../rep-'+str(j+1)+'/Geometry', 'w') as file:
            file.write(opt_coords)
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(natoms):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')
        
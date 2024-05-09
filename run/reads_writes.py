from chemformula import ChemFormula
import numpy as np
def read_xyz(Geometry, atoms, coords):
    with open(Geometry,"r") as file:
        file.readline()
        file.readline()
        for i in range(atoms):
            coords.append(file.readline().split())

def run_script_write(geom_folder):
    with open(geom_folder+"/geom.sh", "w") as file:
        file.write("#$ -cwd -V \n")
        file.write("#$ -l h_vmem=1G,h_rt=1:00:00\n")
        file.write("#$ -pe smp 8\n")
        file.write("module load test qchem\n")
        file.write("module load intel\n")
        file.write("qchem -save -nt 8 opt_freq.inp f.out wf")

def mass_write(EXDIR, atoms, coords):
    with open(EXDIR+"/m.txt", "w") as file:
        for i in range(atoms):
            file.write(str(int(ChemFormula(coords[i][0]).formula_weight))+"\n")

def write_opt_freq_input(geom_folder, atoms, coords):

    with open(geom_folder+"/opt_freq.inp", "w") as file:
        file.write("$molecule\n")
        file.write("0 1\n")
        for i in range(atoms):
            file.write(coords[i][0]+" "+coords[i][1]+" "+coords[i][2]+" "+coords[i][3]+"\n")
        file.write("$end\n")
        file.write("\n")
        file.write("$rem\n")
        file.write("Jobtype		OPT\n")
        file.write("EXCHANGE            BHHLYP !50% HF +  50% Becke88 exchange\n")
        file.write("BASIS               6-31+G*\n")
        file.write("UNRESTRICTED        True\n")
        file.write("MAX_SCF_CYCLES      500\n")
        file.write("SYM_IGNORE          True\n")
        file.write("SCF_Algorithm       DIIS\n")
        file.write("$end\n")
        file.write("\n")
        file.write("@@@\n")
        file.write("\n")
        file.write("$molecule\n")
        file.write("read\n")
        file.write("$end\n")
        file.write("$rem\n")
        file.write("\n")
        file.write("Jobtype	FREQ\n")
        file.write("EXCHANGE            BHHLYP\n")
        file.write("BASIS               6-31+G*\n")
        file.write("$end")

def read_f_out(geom_folder,atoms):
    opt_geoms,modes,mode_cnt=[],[],0
    with open(geom_folder+"/f.out", "r") as file:
        for line in file:
            if line.strip() == "END OF GEOMETRY OPTIMIZER USING LIBOPT3":
                break
        for line in file:
            if line.strip() == "Standard Nuclear Orientation (Angstroms)":
                break
        file.readline()
        file.readline()
        for i in range(atoms):
            num, atom, x, y, z = file.readline().split()
            opt_geoms.append((atom, float(x), float(y), float(z)))
        for line in file:
            if line.strip() =="**        VIBRATIONAL FREQUENCIES (CM**-1) AND NORMAL MODES         **":
                break
        
        while True:
            line=file.readline().split()
            if line==[]:
                continue
            elif line[0]=="Mode:":
                mode_cnt=int(line[-1])
            elif line[0] == "Raman":
                file.readline()
                for i in range(atoms):
                    modes.append(file.readline().split())
            elif line[0] == "*":
                break
    return opt_geoms, modes, mode_cnt

def read_write_breaks(geom_folder,folder,opt_geoms,atoms):
    z_matrix = []
    with open(geom_folder+"/f.out", "r") as file:
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
    outputs = []
    with open(folder+"/bondarr.txt", "w") as file:
        for i in range(atoms):
            for j in range(atoms):
                if bonds[j,i] == 1:
                    file.write((f'{i+1}'+'-'+f'{j+1}'+':'+z_matrix[i][0]+'-'+z_matrix[j][0]+'\n'))

def read_masses(EXDIR):
    m = np.loadtxt(EXDIR+'/m.txt')
    m = m * 1822.887  # Multiply m by 1836
    return m

def write_momentas(opt_geoms,folder,Px,Py,Pz,n,mom_num):
    # Save each momentum to separate files
     for j in range(mom_num):
        with open(folder+f'/Geometry.{j + 1}', 'w') as file:
            for atom in range(n):
               file.write(f'{opt_geoms[atom][0]}  {opt_geoms[atom][1]}  {opt_geoms[atom][2]}  {opt_geoms[atom][3]}\n') 
            file.write("momentum\n")
            # Write Px, Py, and Pz for each atom on the same line
            for atom in range(n):
                # Access the Px, Py, and Pz values using the corresponding indices
                px_value = Px[atom, j]
                py_value = Py[atom, j]
                pz_value = Pz[atom, j]
                file.write(f'{px_value}  {py_value}  {pz_value}\n')


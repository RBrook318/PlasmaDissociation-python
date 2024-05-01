import sys
import os
import subprocess
import numpy as np
from init import Molecule
import init
import qchem as qc
import prop
import output as out

if __name__ == "__main__":
    if len(sys.argv) != 10:
        print("Usage: python script_name.py <reps> <noofcpus> <natoms> <nstates>")
        sys.exit(1)

    try:
        reps= int(sys.argv[1])
        ncpu = int(sys.argv[2])
        n= int(sys.argv[3])
        nstates = int(sys.argv[4])
        nbranch = int(sys.argv[5])
        increment = int(sys.argv[6])
        endstep = int(sys.argv[7])
        restart = str(sys.argv[8])
        geom_start = int(sys.argv[9])
    except ValueError:
        print("Invalid number of CPUs. Please provide a valid integer.")
        sys.exit(1)

if(restart == 'NO'):    
    molecule1 = init.create_molecule(reps+geom_start-1, n,nstates)

    molecule2 = init.create_molecule(None, n,nstates)

    # qc.run_qchem(ncpu,'f.in', molecule1,n, nstates,Guess=False)
    qc.readqchem('f.out',molecule1,9,2)

    out.output_molecule(molecule1)
    startstep = 1

    Guess = True

elif(restart == 'YES'):

    molecule2 = init.create_molecule(None,n,nstates)

    filename = '../output/molecule.json'
    if os.path.exists(filename):

        molecule1 = Molecule.from_json(filename)

        startstep = molecule1.timestep / increment
       
        Guess = False
        qc.run_qchem(ncpu,'f.in', molecule1,n,nstates, Guess=Guess)
        
    else:
        molecule1 = init.create_molecule(reps, n,nstates)
        
        molecule2 = init.create_molecule(None, n,nstates)

        qc.run_qchem(ncpu,'f.in', molecule1,n,nstates, Guess=False)

        out.output_molecule(molecule1)

        startstep = 1

        Guess = True

for i in range(int(startstep), endstep+1):

    molecule2 = prop.prop_1(molecule1, molecule2, n, nstates, increment)
    
    
    
    # qc.run_qchem(ncpu,'f.in', molecule2,n,nstates, Guess=Guess)
    qc.readqchem('f2.out',molecule2,9,2)
    # qc.create_qchem_input('f.in', molecule2, scf_algorithm="DIIS", Guess=True)
    

    molecule1 = prop.prop_2(molecule1, molecule2, n, nstates, increment)

    molecule1, dissociated = prop.fragements(molecule1)
   
    molecule1 = prop.prop_diss(molecule1,increment)

    out.output_molecule(molecule1)

    exit()    
    if dissociated == 0:
        Guess = True
    else:
        Guess = False
    


out.detect_dissociation('../output/xyz.all', 'dissociation.out')
import sys
import os
from init import Molecule
import init
import qchem as qc
import prop
import output as out
import result
import json

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py <reps> <noofcpus> <natoms> <nstates>")
        sys.exit(1)

    if(int(sys.argv[1])==0):
        restart='NO'
    elif(int(sys.argv[1])==1):
        restart='YES'

    with open('../../inputs.json') as f:
        inputs=json.load(f)
    # Check basic arguements
    reps=inputs["setup"]["repeats"]
    ncpu=inputs["setup"]["cores"]
    n=inputs["run"]["Atoms"]
    nstates=inputs["run"]["States"]
    nbranch=inputs["run"]["Branches"]
    increment=inputs["run"]["Timestep"]
    endstep=inputs["run"]["Tot_timesteps"]
    geom_start=inputs["run"]["Geom_start"]
    spin_flip=inputs["run"]["Spin_flip"]

  
if(restart == 'NO'):    
    molecule1 = init.create_molecule(reps+geom_start-1, n,nstates,spin_flip)
    molecule2 = init.create_molecule(None, n,nstates,spin_flip)
    qc.run_qchem(ncpu,'f.in', molecule1,n, nstates,spin_flip,Guess=False)
    out.output_molecule(molecule1)
    startstep = 1
    Guess = True
elif(restart == 'YES'):
    molecule2 = init.create_molecule(None,n,nstates,spin_flip)
    filename = '../output/molecule.json'
    if os.path.exists(filename):
        molecule1 = Molecule.from_json(filename)
        startstep = molecule1.timestep / increment
        Guess = False
        qc.run_qchem(ncpu,'f.in', molecule1,n,nstates,spin_flip, Guess=Guess)
    else:
        molecule1 = init.create_molecule(reps, n,nstates,spin_flip)
        molecule2 = init.create_molecule(None, n,nstates,spin_flip)
        qc.run_qchem(ncpu,'f.in', molecule1,n,nstates,spin_flip, Guess=False)
        out.output_molecule(molecule1)
        startstep = 1
        Guess = True
for i in range(int(startstep), endstep+1):
    molecule2 = prop.prop_1(molecule1, molecule2, n, nstates, increment)
    qc.run_qchem(ncpu,'f.in', molecule2,n,nstates,spin_flip, Guess=Guess)
    molecule1 = prop.prop_2(molecule1, molecule2, n, nstates, increment)
    molecule1, dissociated = prop.fragements(molecule1,spin_flip)
    molecule1 = prop.prop_diss(molecule1,increment)
    out.output_molecule(molecule1)
    if dissociated == 0:
        Guess = True
    else:
        Guess = False
    
result.process_results()
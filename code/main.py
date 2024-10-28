import sys
import os
from init import Molecule
import init
import elec
import prop
import output as out
import result
import json
import time 


start_time = time.time()
if __name__ == "__main__":
    # Load inputs
    with open('../inputs.json') as f:
        inputs=json.load(f)
    # Check basic arguements
    reps=inputs["setup"]["repeats"]
    ncpu=inputs["setup"]["cores"]
    nstates=inputs["run"]["States"]
    nbranch=inputs["run"]["Branches"]
    mult = inputs["run"]["Multiplicity"]
    increment=inputs["run"]["Timestep"]
    endstep=inputs["run"]["Tot_timesteps"]
    geom_start=inputs["run"]["Geom_start"]
    spin_flip=inputs["run"]["Spin_flip"]
    method = inputs["run"]["method"]
    # Check if run is a restart
    if os.path.exists("output/xyz.all"):
        with open("output/xyz.all", "r") as xyz_file:
            lines = xyz_file.readlines()
        for line in reversed(lines):
            if "Timestep:" in line:
                break
        line=line.split()
        third_last_line = int(line[1])
        third_last_line=third_last_line/increment
        print("Time step: ", third_last_line)
        if(third_last_line>=endstep):
            sys.exit("Run already completed")
        else:
            restart='YES'
    else: 
         restart='NO'

if(restart == 'NO'):    
    molecule1 = init.initialize_structure(nstates,spin_flip,mult)
    n = len(molecule1.symbols)
    molecule2 = init.create_empty_molecule(n,nstates,spin_flip)
    molecule1 = elec.run_elec_structure(molecule1, ncpu,n,nstates,spin_flip,method, Guess=False)
    out.output_molecule(molecule1)
    startstep = 1
    Guess = True
elif(restart == 'YES'):
    filename = 'output/molecule.json'
    if os.path.exists(filename):
        molecule1 = Molecule.from_json(filename)
        n = len(molecule1.symbols)
        molecule2 = init.create_empty_molecule(n,nstates,spin_flip)
        startstep = molecule1.timestep / increment
        Guess = False
        molecule1 =elec.run_elec_structure(molecule1, ncpu,n,nstates,spin_flip, method,Guess=False)
    else:
        molecule1 = init.initialize_structure(nstates,spin_flip,mult)
        n = len(molecule1.symbols)
        molecule2 = init.create_empty_molecule(n,nstates,spin_flip)
        molecule1 = elec.run_elec_structure(molecule1, ncpu,n,nstates,spin_flip, method,Guess=False)
        out.output_molecule(molecule1)
        startstep = 1
        Guess = True
for i in range(int(startstep), endstep+1):
    molecule2 = prop.prop_1(molecule1, molecule2, n, nstates, increment)
    molecule2 = elec.run_elec_structure(molecule2, ncpu,n,nstates,spin_flip,method,Guess=Guess)
    molecule1.elecinfo = molecule2.elecinfo
    
    molecule1 = prop.prop_2(molecule1, molecule2, n, nstates, increment)
    molecule1, dissociated = prop.fragements(molecule1,spin_flip)
    molecule1 = prop.prop_diss(molecule1,increment)
    out.output_molecule(molecule1)
    if dissociated == 0:
        Guess = True
    else:
        Guess = False
end_time=time.time()  
with open("output/time.out", "w") as time:
    time.write(f"{(end_time-start_time)/3600}\n")
    time.write(f"{(end_time-start_time)/60}\n")
    time.write(f"{(end_time-start_time)}\n")
result.process_results()


#########################################################################################
#
#   Python Run script for the submission of dissociation after electron impact
#   Written by R Brook                                                 26.07.23
#   Inspired heavily by Oliver Bramley's run script for MCE/CCS
# 
#   Written to make submission of jobs to the HPC easier and quicker
#     
#
#########################################################################################
import sys
import socket
import os
import subprocess
import getpass
import shutil
import json

#########################################################################################
#                * NO NEED TO SCROLL FURTHER IF USING AS BLACKBOX *                     #
#########################################################################################


if __name__=="__main__":
    with open('inputs.json') as f:
        inputs=json.load(f)

    #Check basic arguements
    if(isinstance(inputs["HPC"]["repeats"],int)==False):
        sys.exit("Number of repeats must be an integer")
    elif(isinstance(inputs["HPC"]["cores"],int)==False):
        sys.exit("Number of parallel cores must be an integer")
    elif(inputs["HPC"]["repeats"]<1):
        sys.exit("Not enough runs selected. Must be 1 or greater")
    elif(inputs["HPC"]["cores"]>16):
        sys.exit("Too many cores selected. Maximum of 8 available")
    elif(inputs["HPC"]["cores"]<1):
        sys.exit("Not enough cores selected. Must be 1 or greater")

 
    print("Arguments checked")
    # Get the hostname
    hostname = socket.gethostname()

    # Initialize HPCFLG
    HPCFLG = 0

    # Check for specific substrings in the hostname
    if "arc4" in hostname:
        HPCFLG = "arc"
    elif "arc3" in hostname:
        HPCFLG = "arc"
    elif "aire" in hostname:
        HPCFLG = "aire"

    print(f"HPCFLG: {HPCFLG}")

    if(HPCFLG==0):
        EXDIR="../EXEC"
        if not os.path.exists("../EXEC"):
            os.mkdir("../EXEC")  
    elif(HPCFLG=="arc"):
        EXDIR="/nobackup/"+getpass.getuser()
        HOME_DIR = os.getcwd()
        if not os.path.exists(EXDIR):
            os.mkdir(EXDIR)
    elif(HPCFLG=="aire"):
        EXDIR="/mnt/scratch/"+getpass.getuser()
        if not os.path.exists(EXDIR):
            os.mkdir(EXDIR)
    
    if os.path.exists(EXDIR+"/"+inputs["HPC"]["Runfolder"]):
        value=input("File already exists do you want to delete it? y/n\n")
        if(value=='y'):
            shutil.rmtree(EXDIR+"/"+inputs["HPC"]["Runfolder"])
        else:
            sys.exit("Runfolder already exists. Change the Runfolder name or delte/move it")
    
    os.mkdir(EXDIR+"/"+inputs["HPC"]["Runfolder"])
    EXDIR1=EXDIR+"/"+inputs["HPC"]["Runfolder"]  
    os.mkdir(EXDIR1+"/results")
    os.mkdir(EXDIR1+"/repetitions")
    os.mkdir(EXDIR1+"/results/bonds")
    with open(EXDIR1+"/results/bonds/allbonds.out", 'w'):
        pass
    os.mkdir(EXDIR1+"/results/specifics")
    os.mkdir(EXDIR1+"/results/graphs")
    os.mkdir(EXDIR1+"/setup")
    if (inputs["Molecule_data"]["Geom_flg"]) == "Initial":
        shutil.copy2("../Initial_Geometry/"+inputs["Molecule_data"]["Molecule"]+".txt",EXDIR1+"/setup/initial_guess.txt")
    os.mkdir(EXDIR1+"/setup/tmp")
    shutil.copy2("restart.py",EXDIR1)
    shutil.copy2("inputs.json",EXDIR1)

    for i in range(inputs["HPC"]["repeats"]):
        os.mkdir(EXDIR1+"/repetitions/rep-"+str(i+1))
        os.mkdir(EXDIR1+"/repetitions/rep-"+str(i+1)+"/output")
        os.mkdir(EXDIR1+"/repetitions/rep-"+str(i+1)+"/tmp")
        os.mkdir(EXDIR1+"/repetitions/rep-"+str(i+1)+"/checks")
        if(inputs["Molecule_data"]["Geom_flg"] == "Full"):
            shutil.copy2("../"+inputs["Molecule_data"]["Molecule"]+"/Geom/Geometry."+str(i+inputs["Molecule_data"]["Geom_file_start"]),EXDIR1+"/rep-"+str(i+1)+"/Geometry")
    
    if HPCFLG==1:
            shutil.copytree("../code", EXDIR1+'/code')
    if(inputs["Molecule_data"]["Geom_flg"] ==0):
        shutil.copy2("../"+inputs["Molecule_data"]["Molecule"]+"/bondarr.txt",EXDIR1+"/results")
    
    os.chdir(EXDIR1)
    EXDIR1=os.getcwd()
    if HPCFLG == "arc":
        # Makes Job Submission sctipt for initial setup
        if inputs["Molecule_data"]["Geom_flg"] in ["PubChem", "Initial"]:
            file2="Setup_"+inputs["HPC"]["Runfolder"]+".sh"
            f=open(file2,"w")
            f.write("#$ -cwd -V \n")
            f.write("#$ -l h_vmem=1G,h_rt=02:00:00 \n")
            f.write("#$ -N Setup_"+inputs["HPC"]["Runfolder"]+" \n")
            f.write("#$ -pe smp "+str(inputs["HPC"]["cores"])+" \n") #Use shared memory parallel environemnt 
            if(inputs["Elec_structure"]["method"]=="QChem"):
                f.write("module load  test qchem \n")
                f.write("mkdir $TMPDIR/qchemlocal\n")
                f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
                f.write('qchemlocal=$TMPDIR/qchemlocal\n')
                f.write('export QCHEM_HOME="$qchemlocal"\n')
                f.write('export QC="$qchemlocal"\n')
                f.write('export QCAUX="$QC/qcaux"\n')
                f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
                f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
                f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
                f.write("export QCSCRATCH="+EXDIR1+"/setup/tmp \n")
            f.write("cd "+EXDIR1+"/setup \n")
            f.write("python "+HOME_DIR+"/../code/setup.py")
            f.close()
            command = ['qsub','-N','Setup_'+inputs["HPC"]["Runfolder"], file2]
            subprocess.call(command)
        
        
        file1="Plasma_"+inputs["HPC"]["Runfolder"]+"_1.sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_vmem=8G,h_rt="+str(inputs["HPC"]["hours"])+":00:00 \n")
        f.write("#$ -N Plasma_"+inputs["HPC"]["Runfolder"]+"_1 \n")
        f.write("#$ -pe smp "+str(inputs["HPC"]["cores"])+" \n") #Use shared memory parallel environemnt 
        f.write("#$ -t 1-"+str(inputs["HPC"]["repeats"])+" \n")
        if(inputs["Elec_structure"]["GPU"]==1):
            f.write('#$ -l coproc_p100=1 \n')
        if(inputs["Elec_structure"]["method"]=="QChem"):
            f.write("module load test qchem \n")
            f.write("module load qchem \n")
            f.write("mkdir $TMPDIR/qchemlocal\n")
            f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
            f.write('qchemlocal=$TMPDIR/qchemlocal\n')
            f.write('export QCHEM_HOME="$qchemlocal"\n')
            f.write('export QC="$qchemlocal"\n')
            f.write('export QCAUX="$QC/qcaux"\n')
            f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
            f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
            f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
            f.write("export QCSCRATCH="+EXDIR1+"/repetitions/rep-$SGE_TASK_ID \n")
        f.write("unset GOMP_CPU_AFFINITY KMP_AFFINITY \n")    
        f.write("cd "+EXDIR1+"/repetitions/rep-$SGE_TASK_ID \n")
        f.write("python "+HOME_DIR+"/../code/main.py")
        f.close()
        command = ['qsub','-N','Plasma_'+inputs["HPC"]["Runfolder"]+'_1', '-hold_jid', 'Setup_'+inputs["HPC"]["Runfolder"], file1]
        subprocess.call(command)

        for i in range(inputs["HPC"]["initre"]):
            file1="Plasma_"+inputs["HPC"]["Runfolder"]+"_"+str(i+2)+".sh"
            f=open(file1,"w")
            f.write("#$ -cwd -V \n")
            f.write("#$ -l h_vmem=2G,h_rt=48:00:00 \n")
            f.write("#$ -N Plasma_"+inputs["HPC"]["Runfolder"]+"_"+str(i+2)+" \n")
            f.write("#$ -pe smp "+str(inputs["HPC"]["cores"])+" \n") #Use shared memory parallel environemnt 
            f.write("#$ -t 1-"+str(inputs["HPC"]["repeats"])+" \n")
            f.write("module load qchem \n")
            f.write("mkdir $TMPDIR/qchemlocal\n")
            if(inputs["Elec_structure"]["method"]=="QChem"):
                f.write("mkdir $TMPDIR/qchemlocal\n")
                f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
                f.write('qchemlocal=$TMPDIR/qchemlocal\n')
                f.write('export QCHEM_HOME="$qchemlocal"\n')
                f.write('export QC="$qchemlocal"\n')
                f.write('export QCAUX="$QC/qcaux"\n')
                f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
                f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
                f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
                f.write("export QCSCRATCH="+EXDIR1+"/setup/tmp \n")
            f.write("cd "+EXDIR1+"/repetitions/rep-$SGE_TASK_ID \n")
            f.write("./../main")
            f.close()
            command = ['qsub','-N','Plasma_'+inputs["HPC"]["Runfolder"]+'_'+str(i+2), '-hold_jid', 'Plasma_'+inputs["HPC"]["Runfolder"]+'_'+str(i+1), file1]
            subprocess.call(command)
    elif HPCFLG == "aire":
    # AIRE job submission script
        if inputs["Molecule_data"]["Geom_flg"] in ["PubChem", "Initial"]:
            file2 = "Setup_" + inputs["HPC"]["Runfolder"] + ".sh"
            with open(file2, "w") as f:
                f.write("#!/bin/bash\n")
                f.write("#SBATCH --job-name=Setup_" + inputs["HPC"]["Runfolder"] + "\n")
                f.write("#SBATCH --time=02:00:00\n")
                f.write("#SBATCH --mem=1G\n")
                f.write("#SBATCH --cpus-per-task=" + str(inputs["HPC"]["cores"]) + "\n")
                if inputs["Molecule_data"]["method"] == "QChem":
                    f.write("export LM_LICENSE_FILE=27000@uol-lnx-lic01.leeds.ac.uk\n")
                    f.write("mkdir $TMPDIR/qchemlocal\n")
                    f.write("tar -xzvf /users/" + getpass.getuser() + "/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
                    f.write('qchemlocal=$TMPDIR/qchemlocal\n')
                    f.write('export QCHEM_HOME="$qchemlocal"\n')
                    f.write('export QC="$qchemlocal"\n')
                    f.write('export QCAUX="$QC/qcaux"\n')
                    f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
                    f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
                    f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
                    f.write("export QCSCRATCH=" + EXDIR1 + "/setup/tmp \n")
                f.write("cd " + EXDIR1 + "/setup \n")
                f.write("python /users/"+getpass.getuser()+"/PlasmaDissociation-python/code/setup.py")
            setup_command = ['sbatch', '--parsable', file2]
            setup_job_id = subprocess.check_output(setup_command).strip().decode('utf-8')
            print(setup_job_id)

        file1="Plasma_"+inputs["HPC"]["Runfolder"]+"_1.sh"
        f=open(file1,"w")
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=Plasma_"+inputs["HPC"]["Runfolder"]+"_1 \n")
        f.write("#SBATCH --mem=1G\n")
        f.write("#SBATCH --time="+str(inputs["HPC"]["hours"])+":00:00\n")
        f.write("#SBATCH --cpus-per-task="+str(inputs["HPC"]["cores"])+" \n")
        f.write("#SBATCH --array=1-"+str(inputs["HPC"]["repeats"])+" \n")
        if(inputs["Elec_structure"]["method"]=="QChem"):
            f.write("export LM_LICENSE_FILE=27000@uol-lnx-lic01.leeds.ac.uk\n")
            f.write("mkdir $TMPDIR/qchemlocal\n")
            f.write("tar -xzvf /users/" + getpass.getuser() + "/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
            f.write('qchemlocal=$TMPDIR/qchemlocal\n')
            f.write('export QCHEM_HOME="$qchemlocal"\n')
            f.write('export QC="$qchemlocal"\n')
            f.write('export QCAUX="$QC/qcaux"\n')
            f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
            f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
            f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
            f.write("export QCSCRATCH="+EXDIR1+"/rep-$SLURM_ARRAY_TASK_ID \n")  
        f.write("cd "+EXDIR1+"/repetitions/rep-$SLURM_ARRAY_TASK_ID \n")
        f.write("python /users/"+getpass.getuser()+"/PlasmaDissociation-python/code/main.py")
        f.close()
        plasma_command = ['sbatch', '--dependency=afterok:' + setup_job_id, file1]
        subprocess.call(plasma_command)
    else: 
        for i in range(inputs["HPC"]["repeats"]):
            os.chdir(EXDIR1+"/rep-"+str(i+1))
            if inputs["Molecule_data"]["Geom_flg"] == 1:
                subprocess.call(["python", "../../../code/setup.py"])
            subprocess.call(["python", "../../../code/main.py"])

       
        

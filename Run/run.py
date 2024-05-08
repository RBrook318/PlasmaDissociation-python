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
import random
import shutil
import glob
import csv
import inputs

#########################################################################################
#                              VARIABLES TO SET FOR SIMULATION                          #
#########################################################################################

# Number of repeats 
repeats=5
# #Number of parallel cores per folder/node (max 8)
cores=8
# Name of running folder 
# Default : <method>-<system>-<random number> ie CCS-HP-31254
# Otherwise:  <method>-<system>-<runfolder string>
Runfolder='C2H4-SCF-python'
# Restart Flag 
restart = 'NO'
# Initial restarts 
Initre = 0
#########################################################################################
#                                   END OF INPUTS                                       #
#########################################################################################
#                * NO NEED TO SCROLL FURTHER IF USING AS BLACKBOX *                     #
#########################################################################################


if __name__=="__main__":
    #Check basic arguements
    if(isinstance(repeats,int)==False):
        sys.exit("Number of repeats must be an integer")
    elif(isinstance(cores,int)==False):
        sys.exit("Number of parallel cores must be an integer")
    elif(repeats<1):
        sys.exit("Not enough runs selected. Must be 1 or greater")
    elif(cores>16):
        sys.exit("Too many cores selected. Maximum of 8 available")
    elif(cores<1):
        sys.exit("Not enough cores selected. Must be 1 or greater")

    if(restart=="NO"):
        if(inputs.Geom not in{0,1}):
            sys.exit("Geometry flag must be zero or 1")
        else:
            print("Arguments checked")
            Hostname=socket.gethostname()
            if(Hostname==("login1.arc4.leeds.ac.uk")or(Hostname==("login2.arc4.leeds.ac.uk"))):
                HPCFLG=1
            elif(Hostname==("login1.arc3.leeds.ac.uk")or(Hostname==("login2.arc3.leeds.ac.uk"))):
                HPCFLG=1 
            else:
                HPCFLG=0

        #Makes execution folder and run folder
        if(HPCFLG==0): #change this to 0 before uploading.
            if not os.path.exists("../EXEC"):
                os.mkdir("../EXEC")
            EXDIR="../EXEC"
        else:
            # subprocess.run(['module','load','mkl'])
            os.environ['LOGNAME']
            EXDIR="/nobackup/"+getpass.getuser()

        
        if os.path.exists(EXDIR+"/"+Runfolder):
            value=input("File already exists do you want to delete it? y/n\n")
            if(value=='y'):
                shutil.rmtree(EXDIR+"/"+Runfolder)
            else:
                sys.exit("Runfolder already exists. Change the Runfolder name or delte/move it")
        
        os.mkdir(EXDIR+"/"+Runfolder)
        
        EXDIR1=EXDIR+"/"+Runfolder  
        
        os.mkdir(EXDIR1+"/"+"results")

            
        #Copies input files
        shutil.copy2("run.py",EXDIR1)
        shutil.copy2("inputs.py",EXDIR1)
        shutil.copy2("../"+inputs.Molecule+"/bondarr.txt",EXDIR1+"/"+"results")


        for i in range(repeats):
            path=os.path.join(EXDIR1,"run-"+str(i+1))
            os.mkdir(EXDIR1+"/run-"+str(i+1))
            os.mkdir(EXDIR1+"/run-"+str(i+1)+"/output")
            shutil.copytree("../Code", EXDIR1+"/run-"+str(i+1)+"/Code") 
            shutil.copy2("../"+inputs.Molecule+"/Geom/Geometry."+str(i+inputs.Geom_start),EXDIR1+"/run-"+str(i+1)+"/Code")

                

        os.chdir(EXDIR1)
        EXDIR1=os.getcwd()
        if(HPCFLG==1):
            number=random.randint(99999,1000000)
            file1="Plasma"+str(number)+".sh"
            f=open(file1,"w")
            f.write("#$ -cwd -V \n")
            f.write("#$ -l h_vmem=1G,h_rt=30:00:00 \n")
            f.write("#$ -N file1 \n")
            f.write("#$ -pe smp "+str(cores)+" \n") #Use shared memory parallel environemnt 
            f.write("#$ -t 1-"+str(repeats)+" \n")
            f.write("module load mkl \n")
            f.write("module load test qchem \n")
            f.write("mkdir $TMPDIR/qchemlocal\n")
            f.write('tar -xzvf /nobackup/cm14oab/qchem.tar.gz. -C $TMPDIR/qchemlocal\n')
            f.write('qchemlocal=$TMPDIR/qchemlocal\n')
            f.write('export QCHEM_HOME="$qchemlocal"\n')
            f.write('export QC="$qchemlocal"\n')
            f.write('export QCAUX="$QC/qcaux"\n')
            f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
            f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
            f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
            f.write("module load anaconda \n")
            f.write("source activate base \n")
            f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID/Code \n")
            f.write("python main.py "+"$SGE_TASK_ID"+" "+str(cores)+" "+str(inputs.Atoms)+" "+str(inputs.States)+" "+str(inputs.Branch)+" "+str(inputs.Timestep)+" "+str(inputs.Tot_timesteps)+" "+str(restart)+' '+str(inputs.Geom_start)+' '+str(inputs.Spin_flip))
            f.close()
            # if(cores!=1):
            #     os.environ["OMP_NUM_THREADS"]=str(cores)
            command = ['qsub', '-hold_jid', 'qchemlocal', file1]
            subprocess.call(command)
            # subprocess.call(['qsub', file2])
            for i in range(Initre):
                number=random.randint(99999,1000000)
                file2="Plasma"+str(number)+".sh"
                f=open(file2,"w")
                f.write("#$ -cwd -V \n")
                f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
                f.write("#$ -N file"+str(i+2)+" \n")
                f.write("#$ -pe smp "+str(cores)+" \n") #Use shared memory parallel environemnt 
                f.write("#$ -t 1-"+str(repeats)+" \n")
                f.write("module load mkl \n")
                f.write("module load test qchem \n")
                f.write("module load qchem \n")
                f.write("mkdir $TMPDIR/qchemlocal\n")
                f.write('tar -xzvf /nobackup/cm18rb/qchem.tar.gz. -C $TMPDIR/qchemlocal\n')
                f.write('qchemlocal=$TMPDIR/qchemlocal\n')
                f.write('export QCHEM_HOME="$qchemlocal"\n')
                f.write('export QC="$qchemlocal"\n')
                f.write('export QCAUX="$QC/qcaux"\n')
                f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
                f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
                f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
                f.write("module load anaconda \n")
                f.write("source activate base \n")
                f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID/Code \n")
                f.write("python Main.py "+"$SGE_TASK_ID"+" "+str(cores)+" "+str(inputs.Atoms)+" "+str(inputs.States)+" "+str(inputs.Branch)+" "+str(inputs.Timestep)+" "+str(inputs.Tot_timesteps)+" "+str(restart)+" "+str(inputs.Geom_start)+" "+str(inputs.Spin_flip))
                f.close()
                command = ['qsub', '-hold_jid', 'file'+str(i+1), file2]
                subprocess.call(command)
                
    elif(restart=='YES'):
        print("Arguments checked")
        Hostname=socket.gethostname()
        if(Hostname==("login1.arc4.leeds.ac.uk")or(Hostname==("login2.arc4.leeds.ac.uk"))):
            HPCFLG=1
        elif(Hostname==("login1.arc3.leeds.ac.uk")or(Hostname==("login2.arc3.leeds.ac.uk"))):
            HPCFLG=1 
        else:
            HPCFLG=0
        #If on a SGE machine make job submission file
        EXDIR1=os.getcwd()
        if(HPCFLG==1):
            number=random.randint(99999,1000000)
            file1="Plasma"+str(number)+".sh"
            f=open(file1,"w")
            f.write("#$ -cwd -V \n")
            f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
            f.write("#$ -pe smp "+str(cores)+" \n") #Use shared memory parallel environemnt 
            f.write("#$ -t 1-"+str(repeats)+" \n")
            f.write("module load mkl \n")
            f.write("module load test qchem \n")
            f.write("module load qchem \n")
            f.write("mkdir $TMPDIR/qchemlocal\n")
            f.write('tar -xzvf /nobackup/cm18rb/qchem.tar.gz. -C $TMPDIR/qchemlocal\n')
            f.write('qchemlocal=$TMPDIR/qchemlocal\n')
            f.write('export QCHEM_HOME="$qchemlocal"\n')
            f.write('export QC="$qchemlocal"\n')
            f.write('export QCAUX="$QC/qcaux"\n')
            f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
            f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
            f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
            f.write("module load anaconda \n")
            f.write("source activate base \n")
            f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID/Code \n")
            f.write(" python main.py "+"$SGE_TASK_ID"+" "+str(cores)+" "+str(inputs.Atoms)+" "+str(inputs.States)+" "+str(inputs.Branch)+" "+str(inputs.Timestep)+" "+str(inputs.Tot_timesteps)+" "+str(restart)+' '+str(inputs.Geom_start)+' '+str(inputs.Spin_flip))
            f.close()
            # if(cores!=1):
            #     os.environ["OMP_NUM_THREADS"]=str(cores)
            subprocess.call(['qsub',file1])

        else:
            print('Probably dont run this here')
            

    
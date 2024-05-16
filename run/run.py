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
    if(isinstance(inputs["setup"]["repeats"],int)==False):
        sys.exit("Number of repeats must be an integer")
    elif(isinstance(inputs["setup"]["cores"],int)==False):
        sys.exit("Number of parallel cores must be an integer")
    elif(inputs["setup"]["repeats"]<1):
        sys.exit("Not enough runs selected. Must be 1 or greater")
    elif(inputs["setup"]["cores"]>16):
        sys.exit("Too many cores selected. Maximum of 8 available")
    elif(inputs["setup"]["cores"]<1):
        sys.exit("Not enough cores selected. Must be 1 or greater")

    if(inputs["run"]["Geom_flg"] not in{0,1}):
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


    EXDIR="/nobackup/"+getpass.getuser()+"/scatter"
    
    if os.path.exists(EXDIR+"/"+inputs["setup"]["Runfolder"]):
        value=input("File already exists do you want to delete it? y/n\n")
        if(value=='y'):
            shutil.rmtree(EXDIR+"/"+inputs["setup"]["Runfolder"])
        else:
            sys.exit("Runfolder already exists. Change the Runfolder name or delte/move it")
    
    os.mkdir(EXDIR+"/"+inputs["setup"]["Runfolder"])
    EXDIR1=EXDIR+"/"+inputs["setup"]["Runfolder"]  
    os.mkdir(EXDIR1+"/results")
    os.mkdir(EXDIR1+"/setup")
    os.mkdir(EXDIR1+"/setup/tmp")

    #Copies input files
    shutil.copy2("restart.py",EXDIR1)
    shutil.copy2("inputs.json",EXDIR1)
   
    for i in range(inputs["setup"]["repeats"]):
        os.mkdir(EXDIR1+"/rep-"+str(i+1))
        os.mkdir(EXDIR1+"/rep-"+str(i+1)+"/output")
        os.mkdir(EXDIR1+"/rep-"+str(i+1)+"/tmp")

    os.chdir('../code')
    # subprocess.run(['pyinstaller', 'setup.py', '--onefile'])
    # shutil.copy2("dist/setup",EXDIR1+"/setup/setup")
    shutil.copy2("setup.py",EXDIR1+"/setup/setup.py")
    subprocess.run(['pyinstaller', 'main.py', '--onefile'])
    shutil.copy2("dist/main",EXDIR1+"/main")

    os.chdir(EXDIR1)
    EXDIR1=os.getcwd()
    
    # Makes Job Submission sctipt for iniital setup
    file2="Setup_"+inputs["setup"]["Runfolder"]+".sh"
    f=open(file2,"w")
    f.write("#$ -cwd -V \n")
    f.write("#$ -l h_vmem=1G,h_rt=02:00:00 \n")
    f.write("#$ -N Setup_"+inputs["setup"]["Runfolder"]+" \n")
    f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt 
    f.write("module load qchem \n")
    f.write("mkdir $TMPDIR/qchemlocal\n")
    f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/scatter/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
    f.write('qchemlocal=$TMPDIR/qchemlocal/apps/applications/qchem/6.0.1/1/default\n')
    f.write('export QCHEM_HOME="$qchemlocal"\n')
    f.write('export QC="$qchemlocal"\n')
    f.write('export QCAUX="$QC/qcaux"\n')
    f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
    f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
    f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
    f.write("export QCSCRATCH="+EXDIR1+"/setup/tmp \n")
    f.write("cd "+EXDIR1+"/setup \n")
    f.write("./setup")
    f.close()
    command = ['qsub','-N','Setup_'+inputs["setup"]["Runfolder"], file2]
    subprocess.call(command)
    
   
    file1="Plasma_"+inputs["setup"]["Runfolder"]+"_1.sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
    f.write("#$ -N Plasma_"+inputs["setup"]["Runfolder"]+"_1 \n")
    f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt 
    f.write("#$ -t 1-"+str(inputs["setup"]["repeats"])+" \n")
    f.write("module load qchem \n")
    f.write("mkdir $TMPDIR/qchemlocal\n")
    f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/scatter/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
    f.write('qchemlocal=$TMPDIR/qchemlocal/apps/applications/qchem/6.0.1/1/default\n')
    f.write('export QCHEM_HOME="$qchemlocal"\n')
    f.write('export QC="$qchemlocal"\n')
    f.write('export QCAUX="$QC/qcaux"\n')
    f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
    f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
    f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
    f.write("export QCSCRATCH="+EXDIR1+"/rep-$SGE_TASK_ID/tmp \n")
    f.write("cd "+EXDIR1+"/rep-$SGE_TASK_ID \n")
    f.write("./../main")
    f.close()
    command = ['qsub','-N','Plasma_'+inputs["setup"]["Runfolder"]+'_1', '-hold_jid', 'Setup_'+inputs["setup"]["Runfolder"], file1]
    subprocess.call(command)

    for i in range(inputs["setup"]["initre"]):
        file1="Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i+2)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
        f.write("#$ -N Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i+2)+" \n")
        f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt 
        f.write("#$ -t 1-"+str(inputs["setup"]["repeats"])+" \n")
        f.write("module load qchem \n")
        f.write("mkdir $TMPDIR/qchemlocal\n")
        f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/scatter/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
        f.write('qchemlocal=$TMPDIR/qchemlocal/apps/applications/qchem/6.0.1/1/default\n')
        f.write('export QCHEM_HOME="$qchemlocal"\n')
        f.write('export QC="$qchemlocal"\n')
        f.write('export QCAUX="$QC/qcaux"\n')
        f.write('export QCPROG="$QC/exe/qcprog.exe"\n')
        f.write('export QCPROG_S="$QC/exe/qcprog.exe_s"\n')
        f.write('export PATH="$PATH:$QC/exe:$QC/bin"\n')
        f.write("export QCSCRATCH="+EXDIR1+"/run-$SGE_TASK_ID/tmp \n")
        f.write("cd "+EXDIR1+"/run-$SGE_TASK_ID \n")
        f.write("./../main")
        f.close()
        command = ['qsub','-N','Plasma_'+inputs["setup"]["Runfolder"]+'_'+str(i+2), '-hold_jid', 'Plasma_'+inputs["setup"]["Runfolder"]+'_'+str(i+1), file1]
        subprocess.call(command)

       
        
          
        
    
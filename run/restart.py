#########################################################################################
#
#   Python Restart script for the submission of dissociation after electron impact
#   Written by O Bramley                                                 08.05.24
#   Maininly moving stuff from original run script to here
# 
#   Written to make submission of jobs to the HPC easier and quicker
#     
#
#########################################################################################
import json 
import os
import subprocess
import getpass
if __name__=="__main__":
    with open('inputs.json') as f:
        inputs=json.load(f)
    EXDIR1=os.getcwd()
    
    file1="Plasma_"+inputs["setup"]["Runfolder"]+"_1.sh"
    f=open(file1,"w")
    f.write("#$ -cwd -V \n")
    f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
    f.write("#$ -N Plasma_"+inputs["setup"]["Runfolder"]+"_1 \n")
    f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt 
    f.write("#$ -t "+str(inputs["run"]["Geom_start"])+"-"+str(inputs["setup"]["repeats"])+" \n")
    f.write("module load qchem \n")
    f.write("mkdir $TMPDIR/qchemlocal\n")
    f.write("tar -xzvf /nobackup/"+getpass.getuser()+"/scatter/qchem.tar.gz -C $TMPDIR/qchemlocal\n")
    f.write('qchemlocal=$TMPDIR/qchemlocal/apps/applications/qchem/6.0.1/1/default \n')
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
    command = ['qsub','-N','Plasma_'+inputs["setup"]["Runfolder"]+'_1', file1]
    subprocess.call(command)

    for i in range(inputs["setup"]["initre"]):
        file1="Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i+2)+".sh"
        f=open(file1,"w")
        f.write("#$ -cwd -V \n")
        f.write("#$ -l h_vmem=1G,h_rt=48:00:00 \n")
        f.write("#$ -N Plasma_"+inputs["setup"]["Runfolder"]+"_"+str(i+2)+" \n")
        f.write("#$ -pe smp "+str(inputs["setup"]["cores"])+" \n") #Use shared memory parallel environemnt 
        f.write("#$ -t "+str(inputs["run"]["Geom_start"])+"-"+str(inputs["setup"]["repeats"])+" \n")
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
        f.write("export QCSCRATCH="+EXDIR1+"/setup/tmp \n")
        f.write("cd "+EXDIR1+"/rep-$SGE_TASK_ID \n")
        f.write("./../main")
        f.close()
        command = ['qsub','-N','Plasma_'+inputs["setup"]["Runfolder"]+'_'+str(i+2), '-hold_jid', 'Plasma_'+inputs["setup"]["Runfolder"]+'_'+str(i+1), file1]
        subprocess.call(command)
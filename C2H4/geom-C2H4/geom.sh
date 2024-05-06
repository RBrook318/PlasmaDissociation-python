#$ -cwd -V
#$ -l h_vmem=1G,h_rt=1:00:00
#$ -pe smp 8
module load test qchem
module load intel
qchem -save -nt 8 opt_freq.inp f.out wf
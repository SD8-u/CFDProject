#Example job script used in Leeds ARC4 HPC
#$ -V -cwd

#Request time for execution
#$ -l h_rt=01:00:00

#Request number of cores
#$ -pe smp 40

#Get email at start and end of the job
#$ -m be

#Execute job via containerisation
module load apptainer
apptainer run --cleanenv cfd_apptainer.sif 40 7 10 -5 0.1 0.01 geometry/aneurysm.geo
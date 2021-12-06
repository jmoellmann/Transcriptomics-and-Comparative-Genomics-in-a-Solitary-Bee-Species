#!/bin/bash

#SBATCH -J <job_title>
#SBATCH -A m2_jgu-funcpoll

#SBATCH -p smp (for jobs with c <= 40)
#SBATCH -c <number of cpus per task>
#SBATCH -n <number of tasks>
#SBATCH -t <time in mins>
#SBATCH --mem-per-cpu <mem in MiB>

#SBATCH -o <job_title>.%j.out
#SBATCH -e <job_title>.%j.err

#SBATCH --ramdisk=<mem in GiB>

#PROJECT_DIR=/lustre/project/m2_jgu-funcpoll
#RAMDISK_DIR=/localscratch/${SLURM_JOD_ID}/ramdisk

#rm -f <output_file>
#touch <output_file>

#module purge
#module load <module>

#<code>




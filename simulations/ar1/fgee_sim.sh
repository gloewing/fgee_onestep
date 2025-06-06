#!/bin/bash
#SBATCH --job-name fgee_sim      # Set a name for your job. This is especially useful if you
#SBATCH --partition quick     # Slurm partition to use: quick, norm, 
#SBATCH -c 1        # Number of tasks to run 3
#SBATCH --ntasks-per-core=1	# Do not use hyperthreading (this flag typically used for parallel jobs)
#SBATCH --time 0-0:25       # Wall time limit in D-HH:MM
#SBATCH --mem 15000     # Memory limit for each tasks (in MB) # 1500
#SBATCH -o /home/loewingergc/error/fgee_sim.out  # File to which STDOUT will be written, %j inserts jobid
#SBATCH -e /home/loewingergc/error/fgee_sim.err  # File to which STDERR will be written, %j inserts jobid
#SBATCH --array 1-300
#you can find bash in: /home/loewingergc/bash
module load R/4.4.2

Rscript '/home/loewingergc/fgee/fGEE_sim_ar.R' $1
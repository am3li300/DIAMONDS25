#!/bin/bash
#SBATCH -J edge_add_weight_brain_exp3
#SBATCH --time=01-00:00:00
#SBATCH --mem=128GB
#SBATCH --output=MyJob.%j.%N.out  #saving standard output to file, %j=JOBID,%N=NodeName
#SBATCH --error=MyJob.%j.%N.err   #saving standard error to file, %j=JOBID,%N=NodeName
#SBATCH --mail-type=ALL
#SBATCH -p ccgpu
#SBATCH --mail-user=mert.erden@tufts.edu


./add-remove-edges-tissue.py

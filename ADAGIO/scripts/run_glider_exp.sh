#!/bin/bash

# create lambda array
lambdas=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

# create beta array
betas=(10 20 40 80 160 320 640 1280 2560 5120 10240)

lambdas=(0.1)
betas=(10)

for i in ${lambdas[@]}
    do 
        for j in ${betas[@]}
            do 
                job_name="glider_exp_lam_${i}bet_${j}" 
                
                sbatch -J ${job_name} --time=00-12:00:00 --mem=128GB --output=MyJob.%j.%N.out \
                --error=MyJob.%j.%N.err --mail-type=ALL --mail-user=megan.gelement@tufts.edu \
                glider_experiments.py -t /cluster/tufts/cowenlab/mgelem01/E-MTAB-5214/tissue_files \
                --lamb ${i} --beta ${j} -g /cluster/tufts/cowenlab/mgelem01/tmap/data/alzheimers
            done
    done
#!/bin/bash

#SBATCH -N 1 

#SBATCH -c 12 

#SBATCH --mem-per-cpu=8G 

#SBATCH --time=0-60:00:00 # 60 hours 

#SBATCH --output=my.stdout 

#SBATCH --mail-user=pvs13@case.edu 

#SBATCH --mail-type=ALL 

#SBATCH --job-name="fusion_clustering" 

# Put commands for executing job below this line 


export PYTHONUSERBASE=/home/pvs13/.usr/local/python/3.10.4/lib/python3.12

module load Python

python diffusion_map.py



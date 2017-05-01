#!/bin/bash -l
#SBATCH --nodes=1
#SBATCH --mem=100000
#SBATCH --mail-user=matthew.b.scholz@vanderbilt.edu
#SBATCH --tasks-per-node=12
#SBATCH --time=14-00:00:00
#SBATCH --job-name=Mothur

cd $SLURM_SUBMIT_DIR
echo $PWD
mothur ../mothur/allCommands.m


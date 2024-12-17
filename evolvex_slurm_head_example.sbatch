#!/bin/bash
#SBATCH --job-name=evolvex_head
#SBATCH --account=xxxx
#SBATCH --clusters=xxxx
#SBATCH --partition=xxxx
#SBATCH --output=output_%a.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G 
#SBATCH --time=3-00

module --force purge

module load path/to/cluster/partition
module load Python/xxxx

source /full/path/to/venv/bin/activate

/full/path/to/evolvex_executable /full/path/to/evolvex_config.yaml # To obtain the full path of the evolvex_executable, run "which evolvex".
#!/bin/bash
#SBATCH --job-name=trj_dyn_reg
#SBATCH --partition=gpu_v4_short
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --time=00:10:00
#SBATCH --output=dynamics_reg_%j.out

set -e
source /home/e.lomba/bin/setup_trj
example_dir=/sif/users/e.lomba/Trabajo/IPHS-Mie/trj_analysis_fixed/examples/dynamics_overflow
bash "$example_dir/run.sh" /sif/users/e.lomba/Trabajo/IPHS-Mie/trj_analysis_fixed/bin/trj_analysis \
  "$example_dir/work"

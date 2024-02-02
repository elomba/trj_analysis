#!/bin/bash
#SBATCH -J trj_adp
#SBATCH -e error_%j.err
#SBATCH -o output_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --ntasks-per-node=1
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --constraint=avx2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=adiaz@iqf.csic.es

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
source /opt/intel/intel20/bin/compilervars.sh intel64
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/nv_netcdf/lib64
export PATH=$PATH:$HOME/procluster_netcdf

./trj_analysis.exe input.nml

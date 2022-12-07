#!/bin/sh
#SBATCH -N 5
#SBATCH --ntasks-per-node=48
#SBATCH --time=01:00:00
#SBATCH --job-name=axi-unscaled
#SBATCH --error=job.%J.err_node_48
#SBATCH --output=job.%J.out_node_48
#SBATCH --partition=debug

export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

ulimit -s unlimited
. /opt/ohpc/admin/lmod/8.1.18/init/bash

cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1
mpiexec -n $SLURM_NTASKS hype



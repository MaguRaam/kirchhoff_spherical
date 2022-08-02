#!/bin/sh

#SBATCH -N 25 
#SBATCH --ntasks-per-node=48
#SBATCH --time=12:00:00
#SBATCH --job-name=axi-bubble
#SBATCH --error=job.%J.err_node_48
#SBATCH --output=job.%J.out_node_48
#SBATCH --partition=medium

#export I_MPI_PMI_LIBRARY=/root/rpmbuild/BUILD/slurm-20.11.8/contribs/pmi/.libs/libpmi.so

export OMPI_MCA_btl_openib_allow_ib=1
export OMPI_MCA_btl_openib_if_include="mlx5_0:1"

ulimit -s unlimited
. /opt/ohpc/admin/lmod/lmod/init/sh
##. /opt/ohpc/admin/lmod/8.1.18/init/bash
#module load compiler/intel/2020.4.304
module load openmpi/4.1.2
module load spack/0.17
#source /home-ext/apps/spack/share/spack/setup-env.csh
#spack load openmpi fabrics=auto /nncedca


cd $SLURM_SUBMIT_DIR
export OMP_NUM_THREADS=1
mpiexec -n $SLURM_NTASKS build/bubble


#!/bin/bash
#SBATCH -J wave
#SBATCH -o ./%x.%j.%N.out
#SBATCH -D ./
#SBATCH --get-user-env
#SBATCH --clusters=cm4
#SBATCH --partition=cm4_std
#SBATCH --qos=cm4_std
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=112
#SBATCH --export=NONE
#SBATCH --time=00:15:00
#SBATCH --mail-user=ge59xey@mytum.de
#SBATCH --mail-type=ALL

module load slurm_setup
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module load intel-toolkit/2021.4.0
module unload intel-mpi
module load netcdf-hdf5-all/4.7_hdf5-1.10-gcc12-impi

export NETCDF_LIBRARIES=$(nc-config --libs)
export NETCDF_INCLUDES=$(nc-config --includedir)

export CFLAGS="-std=c11"
export CXXFLAGS="-std=c++20"


mpiexec -np 10 ./build/SWE-MPI-Runner -x 1000 -y 1000

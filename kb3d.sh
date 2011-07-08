#!/bin/ksh
#PBS -N KB3D_FFT_test

#PBS -q rq 

#PBS -l cputim_job=00:05:00
#PBS -l memsz_job=4GB
#PBS -l cpunum_job=1

#PBS -e ./kb3d-%r.err
#PBS -o ./kb3d-%r.out

#PBS -m be -M poghosyan@kit.edu
set -x
# Execute program
cd $PBS_O_WORKDIR
rm kb3d.out
./kb3d

#!/bin/bash

## number of cores
#SBATCH --ntasks=500
#SBATCH --nodes=4

#SBATCH --mail-user=ihepburn@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --job-name=BK50_10nm
#SBATCH --partition=compute

## how much memory per core
#SBATCH --mem-per-cpu=3500m

## maximum time for your simulation, in DAY-HOUR:MINUTE:SECOND
#SBATCH --time=4-0:0:0

module use /apps/unit/DeSchutterU/.modulefiles
module load neuron/2025_03_19_8.2.6

mpirun -n 500 nrngui -mpi figNS.hoc 


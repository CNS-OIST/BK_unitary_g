#!/bin/bash

## number of cores
#SBATCH --ntasks=1
#SBATCH --nodes=1

#SBATCH --job-name=PCcomp
#SBATCH --partition=compute

## how much memory per core
#SBATCH --mem-per-cpu=3500m

## maximum time for your simulation, in DAY-HOUR:MINUTE:SECOND
#SBATCH --time=4-0:0:0

module use /apps/unit/DeSchutterU/.modulefiles
module load neuron/2025_03_19_8.2.6

nrnivmodl mod/

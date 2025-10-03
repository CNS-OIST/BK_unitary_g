# BK_STEPS_NEURON
Models and data for the BK unitary conductance study

The genFigure Python scripts plot the main figures 2-5. 


The STEPS folder contains the Python model script BKmodel_axononly.py from which all data can be generated. 

One example run to run in parallel on 8 cores, with seed 123 and with standard BK and CaP values:
mpirun -n 4 python3 BKmodel_axononly.py 123 1 1


The NEURON folder contains the NEURON model scripts. It also contains example slurm job scripts for cluster usage. 

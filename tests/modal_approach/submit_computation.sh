#!/bin/bash

(
mpirun -np 12 ./runComputation.py > runComputation.log;
mpirun -n 12 SU2_SOL SU2Conf.cfg > SU2_SOL.log;
rm -r flow*.plt;
mkdir Results;
cp surface_flow*.plt history* *.log forces* output* Results/;
zip -r Results.zip Results;
) &

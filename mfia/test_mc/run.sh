#!/usr/bin/bash

# loop several times to relax configurations
for i in `seq 1 5`;
do
   mpirun -n 5 ../MfiaMetrop.exe -i ./ProgramOptionsIn.txt
done

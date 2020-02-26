#!/bin/sh
make
echo using 2 processors
mpirun -np 2 ./b_m.out
mpirun -np 2 ./nb_m.out
mpirun -np 2 ./col_m.out

echo using 4 processors
mpirun -np 4 ./b_m.out
mpirun -np 4 ./nb_m.out
mpirun -np 4 ./col_m.out

echo using 8 processors
mpirun -np 8 ./b_m.out
mpirun -np 8 ./nb_m.out
mpirun -np 8 ./col_m.out
#!/bin/bash

if [ ! -d "Data" ]
then
    mkdir Data
fi

mpirun -np 1 ./a.out

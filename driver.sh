#!/bin/bash

#PBS -q C8
#PBS -l select=1
#PBS -N LBT_test
#PBS -e errLBT.out 
#PBS -o outLBT.out 
#PBS -m abe
#PBS -M yuuka.kanakubo@gmail.com

cd $PBS_O_WORKDIR

time ./LBT parameters.dat ../PYTHIA8/jet_shower_parton.dat lp-posi.dat lp-nega.dat


#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:pnn=4:scivybridge
#PBS -l pmem=2gb
#PBS -A open

cd $PBS_O_WORKDIR
./serial test

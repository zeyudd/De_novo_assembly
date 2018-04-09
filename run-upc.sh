#!/bin/bash
#PBS -l walltime=00:15:00
#PBS -l nodes=1:pnn=4:scivybridge
#PBS -l pmem=2gb
#PBS -A open

# Parameters
P=3				# Number of UPC processes to run
INPUT=test		# Path to your input file

cd $PBS_O_WORKDIR
# Run program
#./serial ${INPUT}
upcrun -np $P -shared-heap=1G ./pgen ${INPUT}

# Sort contigs in both output files to compare
#sort serial.out > serial.sorted
#sort pgen.out > pgen.sorted
# diff -q <file1> <file2> will print nothing if both files are equal
# It will say "Files <file1> and <file2> differ" otherwise.
# Hint: remove -q option to show the lines that differ (not recommended when your output file is large)
#diff -q serial.sorted pgen.sorted
#!/bin/bash
exec > logfile_python.csv
for d in */ ; do
	python3 -W ignore runTriangleCountBenchmark.py ${d:: -4}"adj.mmio"  ${d:: -4}"inc.mmio"	

	#$d"adj.mmio"
       	#mpirun - 4 tc_mpi_cuda.bin $d
	#find .|grep ".mmio$"
done

foldername=$(echo $inputfile)

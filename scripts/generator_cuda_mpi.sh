#!/bin/bash
exec > logfile_cuda_mpi.csv
for d in */ ; do
		#python3 runTriangleCountBenchmark.py ${d:: -4}"adj.mmio"  ${d:: -4}"inc.mmio"	

	#$d"adj.mmio"
       	mpirun -n 4 tc_mpi_cuda.bin $d
	#find .|grep ".mmio$"
done

#foldername=$(echo $inputfile)

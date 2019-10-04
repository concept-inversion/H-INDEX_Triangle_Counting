#!/bin/bash


for d in cd /gpfs/alpine/proj-shared/csc289/trianglecounting/*/ ; do
		#python3 runTriangleCountBenchmark.py ${d:: -4}"adj.mmio"  ${d:: -4}"inc.mmio"	

	#$d"adj.mmio"
       	$d
	#find .|grep ".mmio$"
done



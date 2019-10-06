#!/bin/bash
inputfile=$1
#../../miniTri.exe $inputfile 
awk 'NF{NF-=1};1' < $inputfile > delete.tuple
sed '3s/.*/%&/' delete.tuple > $inputfile".tuple"
rm delete.tuple
./gConvu -c 0 -i $inputfile".tuple" -o $inputfile".edge"
./gConvu -c 3 -v 6302 -i $inputfile".edge" -o csr
./gConvu -c 4  -i csr -o csr
./gConvu -c 5  -i csr -o csr

mv csr.adj_rankbydegree adjacent.bin
mv csr.beg_pos_rankbydegree begin.bin
mv $inputfile".edge" edge
#../../tc-ne-cpu ./

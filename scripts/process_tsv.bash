#!/bin/bash
for d in *.tsv; do
    vertex_coun=$(awk -F " " 'BEGIN{max=0}{for(i=1;i<=NF;i++){if(max<$i)max=$i}}END{printf "%d\n", max+1}' $d)
    echo $vertex_coun
    sed -i '1i %%L1' $d 
    sed -i '2i %L2' $d
    sed -i "3i ${vertex_coun}" $d
done

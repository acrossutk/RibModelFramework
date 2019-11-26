#!/bin/bash
for a1 in 1 2 3 4
do
    filename="S.cerevisiae.S288_a=$a1.beyer.phi.fasta"
    echo "simulating a1=$a1, output filename $filename"
    ./../gilchrist_2009_supl_code/CES_1.0/bin/CES -T tRNA/S.cerevisiae.tRNA.tsv -F REF_FASTA/S.cerevisiae.S288.beyer.phi.fasta 


    
done


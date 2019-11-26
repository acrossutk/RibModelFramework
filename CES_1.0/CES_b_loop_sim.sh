#!/bin/bash
for b in 0.0000025 0.000025 0.00025 0.0025 0.025 0.25
do
    echo "simulating b=$b"
    ./bin/CES -T tRNA_files/S.cerevisiae.tRNA.tsv -F fasta/nse2000.fasta -B $b -O output/vary_b/b=$b


    
done

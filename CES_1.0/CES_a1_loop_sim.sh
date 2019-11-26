#!/bin/bash
for a1 in 1 2 3 4 6 8 10 20 30 40 50 60 70 80 90 100 150 200 250 300 350 400
do
    filename="S.cerevisiae.S288_a=$a1.beyer.phi.fasta"
    echo "simulating a1=$a1, output filename $filename"
    ./bin/CES -T tRNA_files/S.cerevisiae.tRNA.tsv -F fasta/S.cerevisiae.S288.beyer.phi.fasta -A1 $a1 -O output/vary_a1/a1=$a1


    
done


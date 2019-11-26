#shell script to vary pop size for CESpop file
for i in 1.36E7 1.36E6 1.36E5 1.36E4 1.36E3 1.36E2
do
    ./bin/CES -F ./fasta/S.cerevisiae.S288.beyer.phi.fasta ./tRNA_files/S.cerevisiae.tRNA.tsv -Ne $i -o "output/out_Ne=$i"
done

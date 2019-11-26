#shell script for varying protein production rate for CESpop
for i in 0.5 1.0 2.0 5.0 10.0
do
    echo $i
    ./bin/CESpop -F ./fasta/S.cerevisiae.S288.beyer.phi.fasta ./tRNA_files/S.cerevisiae.tRNA.tsv -P $i -o "output/out_Phi=$i"
done

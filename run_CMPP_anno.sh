#!/bin/bash
#run_CMPP_anno.sh -i $your_proteins.fasta -d $/path/to/functional/
#A script to run blastp and hmmscan to annotate C(CAZy)M(MEROPS)P(PHI)P(P450) for fungi.

while getopts i:d: flag
do
    case "${flag}" in
        i) FASTA=${OPTARG};;
        d) PATH_TO_FUNCDB=${OPTARG};;

    esac
done

#echo "$PATH_TO_FUNCDB"

blastp -query $FASTA -db $PATH_TO_FUNCDB/merops/merops_scan.lib -evalue 1e-10 -outfmt 6 -max_target_seqs 5 -num_threads 4 -out $FASTA.merops.btb
blastp -query $FASTA -db $PATH_TO_FUNCDB/PHI/phi-base_current.fas -evalue 1e-10 -outfmt 6 -max_target_seqs 5 -num_threads 4 -out $FASTA.phi.btb
hmmscan -o P450.out --tblout $FASTA.p450.htb --noali --cpu 4 -E 1e-5 $PATH_TO_FUNCDB/P450/P450.hmm.txt $FASTA
hmmscan -o CAZy.out --tblout $FASTA.cazy.htb --noali --cpu 4 -E 1e-5 $PATH_TO_FUNCDB/CAZy/dbCAN-HMMdb-V9.txt $FASTA

grep -v "#" $FASTA.p450.htb|awk -F " " '{print$3}'|sort|uniq > p450.glist
grep -v "#" $FASTA.cazy.htb|awk -F " " '{print$3}'|sort|uniq > cazy.glist
awk '{print$1}' $FASTA.merops.btb |sort|uniq > merops.glist
awk '{print$1}' $FASTA.phi.btb |sort|uniq > phi.glist

mkdir CMPP_out
mv *glist CMPP_out
mv *htb CMPP_out
mv *btb CMPP_out
rm P450.out
rm CAZy.out

echo "Done Blastp and hmmscan!"

python3 CMPP_plot.py $FASTA

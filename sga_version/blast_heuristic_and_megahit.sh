HEURISTICDIR=$1
#LOG=$OUTDIR/$NAME.log
SUF="_filtered_fc2"

#python ../part2_assembly_and_diagnose/why_not_mapping.py $OUTDIR/longest $OUTDIR/not_mapping_longest 21

BLAST=/home/julia/lib/ncbi-blast-2.7.1+/bin/

$BLAST/makeblastdb -in $OUTDIR/megahit/megahit.fa -parse_seqids -dbtype nucl -title megahit

FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
$BLAST/blastn -query ${HEURISTICDIR}/longestfc$SUF.fa -db "$OUTDIR/megahit/megahit.fa" -out ${HEURISTICDIR}/blast_longestfc_to_megahit.tsv     -evalue '1e-5'      -outfmt "$FORMAT"




# megahit to longestfc
python rename_fasta.py ${HEURISTICDIR}/longestfc.fa ${HEURISTICDIR}/longestfc_renamed.fa
$BLAST/makeblastdb -in ${HEURISTICDIR}/longestfc_renamed.fa -parse_seqids -dbtype nucl  -title longestfc

FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
$BLAST/blastn -query "$OUTDIR/megahit/megahit.fa"  -db ${HEURISTICDIR}/longestfc_renamed.fa -out ${HEURISTICDIR}/blast_megahit_to_longestfc.tsv     -evalue '1e-5'      -outfmt "$FORMAT"



BLAST=/home/julia/lib/ncbi-blast-2.7.1+/bin
DIR=$1
SUF="_filtered_fc2"

# download and build database
DBDIR=/mnt/chr5/data/julia/BLASTDB
#mkdir -p $DBDIR
#cd $DBDIR; wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt*.tar.gz
#for file in $DBDIR/nt*tar.gz;
#do
#    cd $DBDIR; tar -xzf $file
#done

FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"


# blast unique sequences from all heuristics
N=1000
#for NAME in ${DIR}/longest ${DIR}/longestfc ${DIR}/bestfc ${DIR}/sga; do
    ### select N random sequences
    #python select_random.py $NAME$SUF.fa ${NAME}_random$N.fa $N
#done
#echo "selected random $N"


#for NAME in ${PREFIX}_max_len ${PREFIX}_max_len_fc ${PREFIX}_best_fc  ${PREFIX}_sga ${PREFIX}_megahit; do
for NAME in longest longestfc bestfc sga; do
    #$BLAST/blastn -query ${DIR}/${NAME}_random$N.fa -db $DBDIR/nt  -out ${DIR}/blast_${NAME}_random$N.tsv  -evalue '1e-5'  -outfmt "$FORMAT" -max_target_seqs 10 -max_hsps 1
    echo $NAME
done



# blast longest sequences from one heuristic and megahit
N=100
python select_longest.py ${DIR}/longestfc$SUF.fa ${DIR}/longestfc_longest$N.fa $N
python select_longest.py ${OUTDIR}/megahit/megahit.fa ${OUTDIR}/megahit/megahit_longest$N.fa $N
echo "selected longest $N"


for NAME in longestfc; do
    FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
    $BLAST/blastn -query ${DIR}/${NAME}_longest$N.fa -db $DBDIR/nt  -out ${DIR}/blast_${NAME}_longest$N.tsv  -evalue '1e-5' -outfmt "$FORMAT" -max_target_seqs 10 -max_hsps 1
    echo $NAME
done
$BLAST/blastn -query ${OUTDIR}/megahit/megahit_longest$N.fa -db $DBDIR/nt  -out  ${OUTDIR}/megahit/megahit_longest$N.tsv  -evalue '1e-5' -outfmt "$FORMAT" -max_target_seqs 10 -max_hsps 1

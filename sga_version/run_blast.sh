
BLAST=/home/julia/lib/ncbi-blast-2.7.1+/bin/blastn
N=1000
PREFIX="my_simplified1_fc2_normed_uniqueM"
#for NAME in ${PREFIX}_max_len ${PREFIX}_max_len_fc ${PREFIX}_best_fc  ${PREFIX}_sga ${PREFIX}_megahit; do
for NAME in ${PREFIX}_max_len_fc ${PREFIX}_megahit; do
    ### select N random sequences
    #python select_random.py $NAME.fa ${NAME}_$N.fa $N
    python select_longest.py $NAME.fa ${NAME}_longest$N.fa $N
done
echo "selected $N"


#for NAME in ${PREFIX}_max_len ${PREFIX}_max_len_fc ${PREFIX}_best_fc  ${PREFIX}_sga ${PREFIX}_megahit; do
for NAME in ${PREFIX}_max_len_fc  ${PREFIX}_megahit; do
    FORMAT="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs"
    #$BLAST -query ${NAME}_$N.fa -db nt -remote        -out blast_$NAME.tsv     -evalue '1e-5'      -outfmt "$FORMAT"    
    $BLAST -query ${NAME}_longest$N.fa -db nt -remote        -out blast_${NAME}_longest$N.tsv     -evalue '1e-5'      -outfmt "$FORMAT"
    echo $NAME
done

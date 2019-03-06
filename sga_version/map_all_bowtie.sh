EXPNAME=$1



SUF="_filtered_fc2"
for NAME in longestfc; do
    GENOMES_FILE=$OUTDIR/heuristics_$EXPNAME/$NAME$SUF.fa
    sh run_bowtie.sh $GENOMES_FILE $OUTDIR/heuristics_$EXPNAME/$NAME
done


GENOMES_FILE=$OUTDIR/sga_scaffold_$EXPNAME/sga-scaffolds_200.fa
NAME="sga"
sh run_bowtie.sh $GENOMES_FILE $OUTDIR/sga_scaffold_$EXPNAME

GENOMES_FILE=$OUTDIR/megahit/megahit_200.fa
NAME="megahit"
sh run_bowtie.sh $GENOMES_FILE $OUTDIR/$NAME

wait






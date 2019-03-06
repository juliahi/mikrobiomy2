EXPNAME=$1


K=21
SUF="_filtered_fc2"

### map to new database with kallisto k

for NAME in longestfc; do
    GENOMES_FILE=$OUTDIR/heuristics_$EXPNAME/$NAME$SUF.fa
    mkdir -p $OUTDIR/heuristics_$EXPNAME/$NAME
    #python analyze_split_fasta.py $GENOMES_FILE
    sh run_kallisto_k_singleend_newkallisto.sh $K $GENOMES_FILE $OUTDIR/heuristics_$EXPNAME/$NAME
    #python why_not_mapping.py $OUTDIR/heuristics_$EXPNAME/$NAME/singleend_newkallisto $OUTDIR/not_mapping_${NAME}_singleend_newkallisto 21
done


python select_contigs.py $OUTDIR/sga_scaffold_$EXPNAME/sga-scaffolds.fa $OUTDIR/sga_scaffold_$EXPNAME/sga-scaffolds_200.fa 200
GENOMES_FILE=$OUTDIR/sga_scaffold_$EXPNAME/sga-scaffolds_200.fa
NAME="sga"
python analyze_split_fasta.py $GENOMES_FILE
sh run_kallisto_k_singleend_newkallisto.sh $K $GENOMES_FILE $OUTDIR/sga_scaffold_$EXPNAME
python why_not_mapping.py $OUTDIR/sga_scaffold_$EXPNAME/singleend_newkallisto $OUTDIR/not_mapping_${NAME}_singleend_newkallisto 21


python select_contigs.py $OUTDIR/megahit/megahit.fa $OUTDIR/megahit/megahit_200.fa 200

GENOMES_FILE=$OUTDIR/megahit/megahit_200.fa
NAME="megahit"
python analyze_split_fasta.py $GENOMES_FILE
sh run_kallisto_k_singleend_newkallisto.sh $K $GENOMES_FILE $OUTDIR/$NAME
#python why_not_mapping.py $OUTDIR/$NAME/singleend_newkallisto $OUTDIR/not_mapping_${NAME}_singleend_newkallisto 21








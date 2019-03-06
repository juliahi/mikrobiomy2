# compile my version of sga
#cd ~/sga_mod/src; ./autogen.sh; ./configure --with-sparsehash=/home/julia/lib/sparsehash/ --with-bamtools=/home/julia/lib/bamtools  --prefix=/home/julia/lib/sga_mod; 
#cd ~/sga_mod/src; make ; make install

SGA_MOD=~/lib/sga_mod/bin/sga

### run sga simplification to compare with my simplification - optional
#name=sga_assemble_simpl_1_200
#$SGA_MOD assemble --cut-terminal=1 --bubble=0 -l 200 -o $name merged.preprocessed_qf5.ec.filter.pass.rmdup_31.asqg >> sga1_200.log
#gunzip -c $name-graph.asqg.gz > $name-graph.asqg


### check if graph is simplified correctly - optional
name=$1
#$SGA_MOD assemble --cut-terminal=0 --bubble=0 -l 200 -o $name  ${name}.asqg >> my_sga1_200.log
#gunzip -c $name-graph.asqg.gz > $name-graph.asqg


SCAFFOLDDIR=$OUTDIR/sga_scaffold_$name
mkdir -p $SCAFFOLDDIR
LOG=$SCAFFOLDDIR/log.log

# Realign reads to the contigs
SGA_NEW=/home/julia/sga_new/src/bin

cp $OUTDIR/${name}_renamed.asqg-contigs.fa $SCAFFOLDDIR/${name}_renamed.asqg-contigs.fa

PRIMARY_CONTIGS=$SCAFFOLDDIR/${name}_renamed.asqg-contigs.fa
PRIMARY_GRAPH=$OUTDIR/${name}_renamed.asqg


BWA_BIN=/home/julia/lib/bwa/bwa

### combine reads in one fasta file
IN1=$OUTDIR/all_r1.fa
IN2=$OUTDIR/all_r2.fa

#touch $IN1
#touch $IN2
#for file in $INDIR/*depl_1.fq.gz; do
#    /home/julia/lib/seqkit fq2fa $file >> $IN1
#done
#for file in $INDIR/*depl_2.fq.gz; do
#    /home/julia/lib/seqkit fq2fa $file >> $IN2
#done


$BWA_BIN index $PRIMARY_CONTIGS 2>> $LOG
$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN1 > $SCAFFOLDDIR/all_r1.sai 2>> $LOG
$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN2 > $SCAFFOLDDIR/all_r2.sai 2>> $LOG
$BWA_BIN sampe $PRIMARY_CONTIGS  $SCAFFOLDDIR/all_r1.sai $SCAFFOLDDIR/all_r2.sai $IN1 $IN2 | samtools view -Sb - > $SCAFFOLDDIR/libPE.bam 2>> $LOG

### Convert the BAM file into a set of contig-contig distance estimates
$SGA_NEW/sga-bam2de.pl -n 1 -m 200 -t 8 --mina 30 -k 100 --prefix $SCAFFOLDDIR/libPE $SCAFFOLDDIR/libPE.bam 2>> $LOG

### Compute copy number estimates of the contigs
$SGA_NEW/sga-astat.py $SCAFFOLDDIR/libPE.bam > $SCAFFOLDDIR/libPE.astat 2>> $LOG

### Build the scaffolds
$SGA_MOD scaffold  -o $SCAFFOLDDIR/scaffolds.scaf -u 1 -c 0.1 -g $PRIMARY_GRAPH --pe $SCAFFOLDDIR/libPE.de $PRIMARY_CONTIGS 2>> $LOG

### Convert the scaffolds to FASTA format#
$SGA_MOD scaffold2fasta --min-gap-length=5 --min-length=200 --use-overlap --distanceFactor=3 --write-unplaced \
    --graph-resolve=best-any --write-names -a $PRIMARY_GRAPH -o $SCAFFOLDDIR/sga-scaffolds.fa $SCAFFOLDDIR/scaffolds.scaf 2>> $LOG



#  best-any: The walk with length closest to the estimated
#                            distance between the contigs will be chosen to resolve the gap.
#                            If multiple best walks are found, the tie is broken arbitrarily.
                     
#  best-unique: as above but if there is a tie no walk will be chosen.

#  unique: only resolve the gap if there is a single walk between the contigs

#  none: do not resolve gaps using the graph

#The most conservative most is unique, then best-unique with best-any being the most aggressive. The default is unique


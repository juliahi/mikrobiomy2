# compile
#cd ~/sga_mod/src; ./autogen.sh; ./configure --with-sparsehash=/home/julia/lib/sparsehash/ --with-bamtools=/home/julia/lib/bamtools  --prefix=/home/julia/lib/sga_mod; 
#cd ~/sga_mod/src; make ; make install

SGA_MOD=~/lib/sga_mod/bin/sga

### run simplification to check - optional
#name=my_full_merged_assemble_simpl_1_200
#$SGA_MOD assemble --cut-terminal=1 --bubble=0 -l 200 -o $name my_full_merged.preprocessed_qf5.ec.filter.pass.rmdup_31.asqg >> sga1_200.log
#gunzip -c $name-graph.asqg.gz > $name-graph.asqg



# run after my simplification to get scaffolds
name=mysimplified1_200_renamed

#$SGA_MOD assemble --cut-terminal=0 --bubble=0 -l 200 -o $name  ${name}.asqg >> my_sga1_200.log



DIR="/mnt/chr7/data/julia/"
PE="/mnt/chr7/data/julia/sga_test_full_notrim_paired_reversed/merged.preprocessed_qf5.ec.filter.pass.rmdup.fa"


SCAFFOLDDIR=$OUTDIR/sga_scaffold
mkdir -p $SCAFFOLDDIR

# Realign reads to the contigs
SGA_NEW=/home/julia/sga_new/src/bin

PRIMARY_CONTIGS=$OUTDIR/$name.asqg-contigs.fa
PRIMARY_GRAPH=$OUTDIR/$name.asqg



IN1=$SCAFFOLDDIR/all_${name}.pe.r1.fa
IN2=$SCAFFOLDDIR/all_${name}.pe.r2.fa
BWA_BIN=/home/julia/lib/bwa/bwa

touch $IN1
touch $IN2

# copy r1 and r2 reads from all files in $PE
for file in $INDIR/*depl_1.fq.gz; do
    /home/julia/lib/seqkit fq2fa $file >> $IN1
done
for file in $INDIR/*depl_2.fq.gz; do
    /home/julia/lib/seqkit fq2fa $file >> $IN2
done


#$BWA_BIN index $PRIMARY_CONTIGS
$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN1 > $IN1.sai
$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN2 > $IN2.sai
$BWA_BIN sampe $PRIMARY_CONTIGS  $IN1.sai $IN2.sai $IN1 $IN2 | samtools view -Sb - > $SCAFFOLDDIR/libPE.bam

### Convert the BAM file into a set of contig-contig distance estimates
$SGA_NEW/sga-bam2de.pl -n 1 -m 200 -t 8 --mina 30 -k 100 --prefix $SCAFFOLDDIR/libPE $SCAFFOLDDIR/libPE.bam

### Compute copy number estimates of the contigs
$SGA_NEW/sga-astat.py $SCAFFOLDDIR/libPE.bam > $SCAFFOLDDIR/libPE.astat

### Build the scaffolds
$SGA_MOD scaffold  -o $SCAFFOLDDIR/scaffolds.scaf -u 1 -c 0.1 -g $PRIMARY_GRAPH --pe $SCAFFOLDDIR/libPE.de $PRIMARY_CONTIGS

### Convert the scaffolds to FASTA format#
$SGA_MOD scaffold2fasta --min-gap-length=5 --min-length=200 --use-overlap --distanceFactor=3 --write-unplaced \
    --graph-resolve=best-any --write-names -a $PRIMARY_GRAPH -o $OUTDIR/sga-scaffolds.fa $SCAFFOLDDIR/scaffolds.scaf



#  best-any: The walk with length closest to the estimated
#                            distance between the contigs will be chosen to resolve the gap.
#                            If multiple best walks are found, the tie is broken arbitrarily.
                     
#  best-unique: as above but if there is a tie no walk will be chosen.

#  unique: only resolve the gap if there is a single walk between the contigs

#  none: do not resolve gaps using the graph

#The most conservative most is unique, then best-unique with best-any being the most aggressive. The default is unique


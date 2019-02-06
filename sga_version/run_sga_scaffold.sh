# compile
#cd ~/sga_mod/src; ./autogen.sh; ./configure --with-sparsehash=/home/julia/lib/sparsehash/ --with-bamtools=/home/julia/lib/bamtools  --prefix=/home/julia/lib/sga_mod; 
#cd ~/sga_mod/src; make ; make install

sga=~/lib/sga_mod/bin/sga

# run simplification to check
name=my_full_merged_assemble_simpl_1_200
#$sga assemble --cut-terminal=1 --bubble=0 -l 200 -o $name my_full_merged.preprocessed_qf5.ec.filter.pass.rmdup_31.asqg >> sga1_200.log 
#gunzip -c $name-graph.asqg.gz > $name-graph.asqg



# run after my simplification to get scaffolds
name=my_simplified1_merged.preprocessed_qf5.ec.filter.pass.rmdup_31_200
#$sga assemble --cut-terminal=0 --bubble=0 -l 200 -o $name  ${name}.asqg >> my_sga1_200.log 
#gunzip -c $name-graph.asqg.gz > $name-graph.asqg

DIR="/mnt/chr7/data/julia/"
PE="/mnt/chr7/data/julia/sga_test_full_notrim_paired_reversed/merged.preprocessed_qf5.ec.filter.pass.rmdup.fa"


# Realign reads to the contigs
SGA_NEW=/home/julia/sga_new/src/bin

#PRIMARY_CONTIGS=${name}-contigs_rename.fa
PRIMARY_CONTIGS=renamed1_200-contigs.fa
PRIMARY_GRAPH=renamed1_200-graph.asqg
IN1=all_${name}.pe.r1.fa
IN2=all_${name}.pe.r2.fa
BWA_BIN=/home/julia/lib/bwa/bwa

#$BWA_BIN index $PRIMARY_CONTIGS
#$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN1 > $IN1.sai
#$BWA_BIN aln -t 16 $PRIMARY_CONTIGS $IN2 > $IN2.sai
#$BWA_BIN sampe $PRIMARY_CONTIGS  $IN1.sai $IN2.sai $IN1 $IN2 | samtools view -Sb - > libPE3.bam

# Convert the BAM file into a set of contig-contig distance estimates
#$SGA_NEW/sga-bam2de.pl -n 1 -m 200 -t 8 --mina 30 -k 100 --prefix libPE3 libPE3.bam
#exit
# Compute copy number estimates of the contigs
#$SGA_NEW/sga-astat.py libPE.bam > libPE.astat

# Build the scaffolds
#$sga scaffold  -o scaffolds3.scaf -u 1 -c 0.1 -g $PRIMARY_GRAPH --pe libPE3.de $PRIMARY_CONTIGS

# Convert the scaffolds to FASTA format#
$sga scaffold2fasta --min-gap-length=5 --min-length=200 --use-overlap --distanceFactor=3 --write-unplaced --graph-resolve=best-any --write-names -a $PRIMARY_GRAPH -o sga-scaffolds4.fa scaffolds2.scaf



#  best-any: The walk with length closest to the estimated
#                            distance between the contigs will be chosen to resolve the gap.
#                            If multiple best walks are found, the tie is broken arbitrarily.
                     
#  best-unique: as above but if there is a tie no walk will be chosen.

#  unique: only resolve the gap if there is a single walk between the contigs

#  none: do not resolve gaps using the graph

#The most conservative most is unique, then best-unique with best-any being the most aggressive. The default is unique


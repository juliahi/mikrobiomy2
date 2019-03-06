### Whole analysis of SGA overlap graph ###

### External software: /home/julia/seqkit, sga, 
### Python: networkx, 

OUTDIR="/home/julia/mikrobiomy_results"
export OUTDIR
mkdir -p $OUTDIR

INDIR="/home/julia/Wyniki_sekwencjonowania"
export INDIR


########################## reverse-complement R2 reads ##############################
#mkdir -p $INDIR/reversed
#for PROBE in '6683_16-06-2015' '6685_04-06-2015' '6685_16-06-2015' '6690_04-06-2015' '6690_16-06-2015' '6695_04-06-2015' #'6695_16-06-2015' '6704_04-06-2015' '6704_16-06-2015'; do
#    echo "reversing $PROBE"
#	#/home/julia/seqkit seq -r -p $INDIR/${PROBE}_depl_2.fq.gz -o $INDIR/reversed/${PROBE}_depl_2.fq.gz &
#done
#wait


######################## Run SGA to get Overlap Graph ###############################
#sh run_sga_separate_full_reverseR2.sh
### results in /mnt/chr7/data/julia/sga_test_full_notrim_paired_reversed
### unchanged graph: merged.preprocessed_qf5.ec.filter.pass.rmdup_31.asqg
### simplified: merged.preprocessed_qf5.ec.filter.pass.rmdup_31-graph.asqg


######################### Analysis of Overlap graph - simplification stages #########
echo "\n#########################################################\n"
NAME=mysimplified1_deadends200_minlen200
echo `date`, "$NAME, running with parameters:"
cat common_parameters.py
echo "-----------------------------------------------------------\n"

#python load_filter_graph.py                   ### graph after renaming in: $OUTDIR/$NAME_renamed.asqg

### visualise stats with jupyter:
## SGA graph analysis - simplification as in SGA


######################### Run heuristics ############################################

### Prepare heuristics and other files (like blast results) for comparison
#sh run_sga_scaffold.sh    $NAME         # output in OUTDIR/sga
#python run_heuristics.py        # output in OUTDIR/heuristics_$NAME/longest longestfc bestfc sga

### Statistics for heuristics
## SGA graph - heuristics - FC=2 normalized

#mkdir -p $OUTDIR/megahit
#cp /mnt/chr4/mikrobiomy-2/megahit_results/all/long_contigs_200.fa $OUTDIR/megahit/megahit.fa
#sh run_blast_local.sh $OUTDIR/heuristics_$NAME

### Compare using mapping with kallisto (paired end and single end) and bowtie2
#sh map_all_kallisto_singleend.sh $NAME
#sh map_all_bowtie.sh $NAME

#python prepare_coverages.py singleend
#python prepare_coverages.py bowtie

##Jupyter - compare heuristics with megahit using kallisto and bowtie mapping
## SGA graph - heuristics - megahit

sh blast_heuristic_and_megahit.sh $OUTDIR/heuristics_$NAME









K=$1

N=$2

#OUTDIR=/mnt/chr4/mikrobiomy-2/small_test_$N/velvet_${K}_track
OUTDIR=/home/julia/dane



INDIR=/mnt/chr4/mikrobiomy-2

mkdir -p $OUTDIR

FQ1=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_04-06-2015_depl_1.fq
FQ2=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_16-06-2015_depl_1.fq

T1=$OUTDIR/small_test_$N/s1.fq
T2=$OUTDIR/small_test_$N/s2.fq

echo $(($N*4)) 
head -n $(($N*4)) $FQ1 > $T1
head -n $(($N*4)) $FQ2 > $T2



OUTDIR=$OUTDIR/small_test_$N/velvet_${K}_track


mkdir -p $OUTDIR
~/velvet_1.2.10/velveth $OUTDIR $K -fastq -short $T1 -short2 $T2 > $OUTDIR/log.txt
wait
~/velvet_1.2.10/velvetg $OUTDIR -cov_cutoff 0 -exp_cov 1  -read_trkg yes  > $OUTDIR/log.txt
#~/velvet_1.2.10/velvetg $OUTDIR -cov_cutoff auto -exp_cov auto  -read_trkg yes  > $OUTDIR/log.txt
wait


# python Graph.py $OUTDIR/LastGraph $OUTDIR/zachlan_ $FQ1,$FQ2 0,1

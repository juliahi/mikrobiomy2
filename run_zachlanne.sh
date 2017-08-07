
INDIR="/mnt/chr4/mikrobiomy-2"
# FQ1=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_04-06-2015_depl_1.fq.gz
# FQ2=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_16-06-2015_depl_1.fq.gz

K=21


#time python zachlan_sumfC.py "$INDIR/small_test/velvet_13_test/LastGraph"  $INDIR/small_test/velvet_13_test/zachlan_  $INDIR/small_test/t.fa 0
#exit

V=0.4
LOG=zachlan_v$V.out
COMMAND="time  python zachlan_sumfC.py"


<<<<<<< HEAD
=======
for N in 1000; # 10000 1000000 ; 
do
>>>>>>> 0036e1502ce33b84f7fe37571996326d4da65453

KALLISTO="/home/julia/kallisto_linux-v0.43.0/kallisto"


for N in 10000 1000000 ; 
do
    DIR=$INDIR/small_test_$N/velvet_${K}_track
# #     now=$(date '+%d/%m/%Y %H:%M:%S');
# #     echo $N, "$now" >>$LOG
# #     echo "velvet"
# #     #time sh small_test.sh $K $N
# #     
# #     echo "Graph"
# #     $COMMAND  -g $DIR/LastGraph  -o $DIR/zachlan_v${V}_ --minlen 200 -f 0,$N 1,$N --profile 2>> $LOG #>> $DIR/zachlan_profile$V.txt
# #     now=$(date '+%d/%m/%Y %H:%M:%S');
# #     echo 'finished', $N, "$now" >> $LOG
# #     
# #     
    ####kallisto k=15 
    OUTDIR=$DIR/velvet_mapping
    mkdir -p $OUTDIR
    FQ=$INDIR/small_test_$N/s1.fq
    K2=15
    INDEX_FILE=$OUTDIR/kallisto_index_${K2}.idx
    GENOMES_FILE=$DIR/contigs.fa
    if [ ! -f  $INDEX_FILE ]
    then
        $KALLISTO index -k ${K2} -i $INDEX_FILE $GENOMES_FILE
    fi
    $KALLISTO quant -i $INDEX_FILE -o $OUTDIR/kallisto${K2} --pseudobam -b 100 --single -l 100 -s 20  $FQ |  samtools view -Sb - > $OUTDIR/pseudoal.bam
    
    ####kallisto na naszych
    OUTDIR=$DIR/zachlan_mapping
    mkdir -p $OUTDIR
    INDEX_FILE=$OUTDIR/kallisto_index_${K2}.idx
    GENOMES_FILE=$DIR/zachlan_v${V}_assemblies.fa
    if [ ! -f  $INDEX_FILE ]
    then
        $KALLISTO index -k ${K2} -i $INDEX_FILE $GENOMES_FILE
    fi
    $KALLISTO quant -i $INDEX_FILE -o $OUTDIR/kallisto${K2} --pseudobam -b 100 --single -l 100 -s 20 $FQ  |  samtools view -Sb - > $OUTDIR/pseudoal.bam
    

    
done

echo "----------------------"


echo "velvet31"
#all data

INDIR="/home/julia/dane"

FILES="1,12066571 0,11557403 1,11887627 0,12888782 1,11275417 0,7374023 1,12136999  0,3018723  1,8376929 "



# # DIR=$INDIR/velvet_31_expcovauto/all
# DIR=$INDIR/velvet_31_expcovauto
#     now=$(date '+%d/%m/%Y %H:%M:%S');
# echo "all $DIR", "$now" >>$LOG
# $COMMAND -g "$DIR/LastGraph"  -o "$DIR/zachlan_v${V}_" -f $FILES --paired --minlen 200 --minfc 4 2>> $LOG #>> $DIR/zachlan_profile$V.txt
#     now=$(date '+%d/%m/%Y %H:%M:%S');
# echo "finished all $DIR", "$now" >>$LOG
# 
# 
# echo "----------------------\n\n"



echo "metavelvet21"
#DIR=$INDIR/metavelvet_21/all
DIR=$INDIR/metavelvet_21
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "all $DIR", "$now" >>$LOG
$COMMAND -g  "$DIR/meta-velvetg.LastGraph"  -o "$DIR/zachlan_v${V}_" -f $FILES --paired --minlen 200 --minfc 4 2>> $LOG #>> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "finished all $DIR", "$now" >>$LOG

wait 

echo "----------------------\n\n"

cat $LOG


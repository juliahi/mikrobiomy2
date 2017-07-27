
INDIR="/mnt/chr4/mikrobiomy-2"
# FQ1=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_04-06-2015_depl_1.fq.gz
# FQ2=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_16-06-2015_depl_1.fq.gz

K=21


#time python zachlan_sumfC.py "$INDIR/small_test/velvet_13_test/LastGraph"  $INDIR/small_test/velvet_13_test/zachlan_  $INDIR/small_test/t.fa 0
#exit

V=0.3
LOG=zachlan_v$V.out
for N in 1000; #10000 100000 1000000 ; 
do

    DIR=$INDIR/small_test_$N/velvet_${K}_track
    now=$(date '+%d/%m/%Y %H:%M:%S');
    echo $N, "$now" >>$LOG
    echo "velvet"
    #time sh small_test.sh $K $N
    
    FQ1=$INDIR/small_test_$N/s1.fq
    FQ2=$INDIR/small_test_$N/s2.fq
    
    
    echo "Graph"
    time  python -m cProfile zachlan_sumfC.py -g $DIR/LastGraph  -o $DIR/zachlan_v0.2_ -f 0,$FQ1 1,$FQ2 2>> $LOG >> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
    echo 'finished', $N, "$now" >> $LOG
    
done

echo "velvet31"
#all data
F=$INDIR/Wyniki_sekwencjonowania/demultiplexed
E=-06-2015_depl_1.fq.gz
###FILES="$F/6683_16$E,$F/6685_04$E,$F/6685_16$E,$F/6690_04$E,$F/6690_16$E,$F/6695_04$E,$F/6695_16$E,$F/6704_04$E,$F/6704_16$E"











FILES="1,12066571 0,11557403 1,11887627 0,12888782 1,11275417 0,7374023 1,12136999  0,3018723  1,8376929 "
DIR=$INDIR/velvet_31_expcovauto/all
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "all $DIR", "$now" >>$LOG
time  python -m cProfile zachlan_sumfC.py -g "$DIR/LastGraph"  -o "$DIR/zachlan_" -f $FILES --paired 2>> $LOG >> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "finished all $DIR", "$now" >>$LOG




echo "metavelvet21"
DIR=$INDIR/metavelvet_21/all
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "all $DIR", "$now" >>$LOG
time  python -m cProfile zachlan_sumfC.py -g  "$DIR/meta-velvetg.LastGraph"  -o "$DIR/zachlan_" -f $FILES --paired 2>> $LOG >> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "finished all $DIR", "$now" >>$LOG

wait 
cat $LOG
#
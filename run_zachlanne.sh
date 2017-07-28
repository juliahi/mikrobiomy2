
INDIR="/mnt/chr4/mikrobiomy-2"
# FQ1=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_04-06-2015_depl_1.fq.gz
# FQ2=$INDIR/Wyniki_sekwencjonowania/demultiplexed/6685_16-06-2015_depl_1.fq.gz

K=21


#time python zachlan_sumfC.py "$INDIR/small_test/velvet_13_test/LastGraph"  $INDIR/small_test/velvet_13_test/zachlan_  $INDIR/small_test/t.fa 0
#exit

V=0.4
LOG=zachlan_v$V.out
COMMAND="time  python zachlan_sumfC.py"


for N in 1000 10000 1000000 ; 
do

    DIR=$INDIR/small_test_$N/velvet_${K}_track
    now=$(date '+%d/%m/%Y %H:%M:%S');
    echo $N, "$now" >>$LOG
    echo "velvet"
    #time sh small_test.sh $K $N
    
    echo "Graph"
    $COMMAND  -g $DIR/LastGraph  -o $DIR/zachlan_v${V}_ --minlen 200 -f 0,$N 1,$N --profile 2>> $LOG #>> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
    echo 'finished', $N, "$now" >> $LOG
    
done

exit


echo "velvet31"
#all data

FILES="1,12066571 0,11557403 1,11887627 0,12888782 1,11275417 0,7374023 1,12136999  0,3018723  1,8376929 "
DIR=$INDIR/velvet_31_expcovauto/all
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "all $DIR", "$now" >>$LOG
$COMMAND -g "$DIR/LastGraph"  -o "$DIR/zachlan_v${V}_" -f $FILES --paired --minlen 200 2>> $LOG #>> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "finished all $DIR", "$now" >>$LOG




echo "metavelvet21"
DIR=$INDIR/metavelvet_21/all
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "all $DIR", "$now" >>$LOG
$COMMAND -g  "$DIR/meta-velvetg.LastGraph"  -o "$DIR/zachlan_v${V}_" -f $FILES --paired --minlen 200 2>> $LOG #>> $DIR/zachlan_profile$V.txt
    now=$(date '+%d/%m/%Y %H:%M:%S');
echo "finished all $DIR", "$now" >>$LOG

wait 
cat $LOG
#
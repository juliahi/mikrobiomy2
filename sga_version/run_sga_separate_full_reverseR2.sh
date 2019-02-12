OUTDIR=/mnt/chr7/data/julia/sga_test_full_notrim_paired_reversed
INDIR=/home/julia/Wyniki_sekwencjonowania

mkdir -p $OUTDIR

LOG=$OUTDIR/log.txt



######################## run SGA ############################
QUALFILTER=5
CORRECT=1
FILTER=0
OVERLAP=31
echo `date` 'start'

logcommand () {
	COMMAND="$2"
	LOGFILE=$1
	echo `date` "running: $COMMAND"  >> $LOGFILE
	$COMMAND  2>>$LOGFILE >>$LOGFILE
}



############# Preprocessing ###############
CORRECTK=21
processprobe () {
    probe=$1
    PREFIX=${probe}
    SUF=""
    SUF="$SUF.preprocessed_qf$QUALFILTER"
    logcommand $OUTDIR/${probe}.log "sga preprocess -o ${OUTDIR}/${PREFIX}$SUF.fa -m 50 --suffix=:$probe --pe-mode=1 --quality-filter=$QUALFILTER --discard-quality $INDIR/${probe}_depl_1.fq.gz $INDIR/reversed/${probe}_depl_2.fq.gz "
    if [ $CORRECT -eq 1 ]; then
	logcommand $OUTDIR/${probe}.log  "sga index -a ropebwt -t 16 --prefix=$OUTDIR/$PREFIX$SUF --no-reverse  $OUTDIR/$PREFIX$SUF.fa"
	echo `date` 'index' $probe
	logcommand $OUTDIR/${probe}.log  "sga correct -t 16 -k $CORRECTK --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa  "
        echo `date` 'correct' $probe
        SUF=${SUF}.ec
	mv $PREFIX$SUF.fa $OUTDIR/
    fi 
    logcommand $OUTDIR/${probe}.log  "sga index -a ropebwt -t 16 --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa"
    echo `date` 'index2' $probe

    if [ $FILTER -eq 1 ]; then
        ######sga filter --no-duplicate-check --no-kmer-check -t 16 $OUTDIR/$PREFIX$SUF.fa 
	logcommand $OUTDIR/${probe}.log  "sga filter --substring-only --kmer-size=$CORRECTK --kmer-threshold=3  -t 16 $OUTDIR/$PREFIX$SUF.fa "
	SUF=${SUF}.filter.pass
	echo `date`  'filter' $probe
    fi

}



for probe in '6683_16-06-2015' '6685_04-06-2015' '6685_16-06-2015' '6690_04-06-2015' '6690_16-06-2015' '6695_04-06-2015' '6695_16-06-2015' '6704_04-06-2015' '6704_16-06-2015'; do
    processprobe $probe &
    sleep 1
done

SUF=".preprocessed_qf$QUALFILTER"
if [ $CORRECT -eq 1 ]; then
        SUF=${SUF}.ec
fi
if [ $FILTER -eq 1 ]; then
        SUF=${SUF}.filter.pass
fi


wait


########### Merging control and treated, rmdup... ########## 
# when more files add merging and output names

CONTR='control'
TREAT='treated'
echo '' > $OUTDIR/$CONTR$SUF.fa
echo '' > $OUTDIR/$TREAT$SUF.fa
for probe in '6685_04-06-2015' '6690_04-06-2015' '6695_04-06-2015'  '6704_04-06-2015' ; do
    cat $OUTDIR/$probe$SUF.fa >> $OUTDIR/$CONTR$SUF.fa
done
for probe in '6683_16-06-2015' '6685_16-06-2015' '6690_16-06-2015' '6695_16-06-2015' '6704_16-06-2015'; do
    cat $OUTDIR/$probe$SUF.fa >> $OUTDIR/$TREAT$SUF.fa
done


#remove duplicates
rmdup () {
	probe=$1
        SUF=$2
	logcommand $OUTDIR/${probe}$SUF.log  "sga rmdup -t 16 --prefix=$OUTDIR/$probe$SUF   $OUTDIR/$probe$SUF.fa "
	mv  ${probe}$SUF.rmdup*  "$OUTDIR/"  		#### bo opcja --prefix nie działa
}

for PREFIX in $CONTR $TREAT; do
	logcommand $OUTDIR/${probe}$SUF.log   "sga index -a ropebwt -t 16  --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa"
	echo `date` "$PREFIX index"
	rmdup $PREFIX $SUF
	echo `date` "$PREFIX rmdup"
done

DUPSUF=$SUF.rmdup

########### Merging treated and control, remove duplicates ##########
PREFIX=merged
#####sga merge usuwa zawarte sekwencje. Dla wersji z trimingiem odczytow trzeba je przygotować inaczej (patrz poniżej -->). Dla odczytów równej długości to jest chyba szybsze.
logcommand $OUTDIR/merged$SUF.log  "sga merge -p $OUTDIR/$PREFIX$SUF -t 8 $OUTDIR/${CONTR}${DUPSUF}.fa $OUTDIR/${TREAT}${DUPSUF}.fa "
echo `date` 'merge'


########## --> wersja dla odczytów różnej długości ?
#PREFIX=catmerged
#cat $OUTDIR/${CONTR}${DUPSUF}.fa > $OUTDIR/$PREFIX$SUF.fa
#cat $OUTDIR/${TREAT}${DUPSUF}.fa >> $OUTDIR/$PREFIX$SUF.fa

###sga index -a ropebwt -t 16 --prefix=$OUTDIR/$PREFIX$SUF --no-reverse  $OUTDIR/$PREFIX$SUF.fa   #rozważyć korekcję tutaj a nie wcześniej
###sga correct -t 16 -k $CORRECTK  --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa    
###SUF=$SUF.correct$CORRECTK
#sga index -a ropebwt -t 16  --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa  
########## end: wersja dla odczytów różnej długości

logcommand $OUTDIR/merged$SUF.log  "sga rmdup -t 8 --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa " 
echo `date` 'rmdup merged' >$LOG  2>$LOG 
SUF=$SUF.rmdup
mv  $PREFIX$SUF*  $OUTDIR/  		#### bo opcja --prefix nie działa


#### finding which duplicate duplicates which sequence using overlap graph build
logcommand $OUTDIR/merged$SUF.log  "sga overlap -m $OVERLAP -t 16 --target-file=$OUTDIR/${PREFIX}$SUF.fa $OUTDIR/${PREFIX}$SUF.dups.fa "
echo `date` 'overlap duplicates' 
mv ${PREFIX}$SUF* $OUTDIR/				#### bo opcja --prefix nie działa
DUPSGRAPH=$OUTDIR/${PREFIX}$SUF.dups.${PREFIX}$SUF.asqg
gunzip -c ${DUPSGRAPH}.gz > $DUPSGRAPH



########## Overlaps ###########

logcommand $OUTDIR/merged$SUF.log  "sga overlap -m $OVERLAP -t 16 --prefix=$OUTDIR/${PREFIX}$SUF $OUTDIR/${PREFIX}$SUF.fa "
echo `date` 'overlap'
mv  ${PREFIX}$SUF.asqg.gz  $OUTDIR/${PREFIX}${SUF}_${OVERLAP}.asqg.gz   		#### bo opcja --prefix nie działa
SUF=${SUF}_${OVERLAP}
gunzip -c $OUTDIR/${PREFIX}$SUF.asqg.gz > $OUTDIR/${PREFIX}$SUF.asqg  



######### Asemble -- get contigs and create simplified chunk graph ##########
logcommand $OUTDIR/merged$SUF.log  "sga assemble -o $OUTDIR/${PREFIX}$SUF $OUTDIR/${PREFIX}$SUF.asqg.gz"
echo `date` 'assemble'
### rozpakowanie "oczyszczonego" grafu nałożeń - assembly wykonuje kilka rzeczy żeby graf uprościć 
gunzip -c $OUTDIR/${PREFIX}$SUF-graph.asqg.gz > $OUTDIR/${PREFIX}$SUF-graph.asqg





#sga scaffold -o $$OUTDIR/${PREFIX}$SUF.scaf --pe=$OUTDIR/${PREFIX}$SUF.fa  $OUTDIR/${PREFIX}$SUF-contigs.fa
#echo `date` 'scaffold'
#sga scaffold2fasta -o $OUTDIR/${PREFIX}$SUF.scaffolds.fa -a $OUTDIR/${PREFIX}$SUF-graph.asqg.gz $OUTDIR/${PREFIX}$SUF.scaf
#echo `date` 'scaffold2fasta'



echo "---------------------------------"
echo "graph:" $OUTDIR/${PREFIX}$SUF.asqg  
echo "duplicates for control:" $OUTDIR/$CONTR$DUPPREF.fa
echo "duplicates for treated:" $OUTDIR/$TREAT$DUPPREF.fa
echo "duplicates for merged:" $DUPSGRAPH

exit



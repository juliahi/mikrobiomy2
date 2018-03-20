L=200 #lenght of contig


OUTDIR=/mnt/chr7/data/julia/sga_test_notrim
INDIR=/home/julia/Wyniki_sekwencjonowania

INFILES1=`echo $INDIR/*depl_?.fq.gz`
echo $INFILES1

mkdir -p $OUTDIR

LOG=$OUTDIR/log.txt

######################## run SGA ############################
QUALFILTER=20
CORRECT=1
FILTER=1
OVERLAP=51
echo `date` 'start'

############# Preprocessing ###############
CORRECTK=21
processprobe () {
    probe=$1
    PREFIX=${probe}
    SUF=""
    #sga preprocess -o ${OUTDIR}/${PREFIX}$SUF.fa -m $OVERLAP --suffix=":$probe" --quality-trim=$QUALFILTER --discard-quality $INDIR/${probe}_depl_1.fq.gz 
    
    SUF="$SUF.preprocessed_q$QUALFILTER"
    if [ $CORRECT -eq 1 ]; then
	#sga index -a ropebwt -t 16 --prefix="$OUTDIR/$PREFIX$SUF" --no-reverse  $OUTDIR/$PREFIX$SUF.fa
	echo `date` 'index' $probe
	#sga correct -t 16 -k $CORRECTK --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa  
        echo `date` 'correct' $probe
        SUF=${SUF}.ec
	mv $PREFIX$SUF.fa $OUTDIR/
    fi 
    #sga index -a ropebwt -t 16 --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa
    echo `date` 'index2' $probe

    if [ $FILTER -eq 1 ]; then
        ######sga filter --no-duplicate-check --no-kmer-check -t 16 $OUTDIR/$PREFIX$SUF.fa 
	#sga filter --substring-only --kmer-size=$CORRECTK --kmer-threshold=3  -t 16 $OUTDIR/$PREFIX$SUF.fa 
	SUF=${SUF}.filter.pass
	echo `date`  'filter' $probe
    fi

}

for probe in '6685_04-06-2015' '6685_16-06-2015'; do
    processprobe $probe &
done
SUF=".preprocessed_q$QUALFILTER"
if [ $CORRECT -eq 1 ]; then
        SUF=${SUF}.ec
fi
if [ $FILTER -eq 1 ]; then
        SUF=${SUF}.filter.pass
fi

wait



########### Merging control and treated, rmdup... ########## 
# when more files add merging and output names
CONTR='6685_04-06-2015'
TREAT='6685_16-06-2015'

#remove duplicates
rmdup () {
	probe=$1
        SUF=$2
	#sga rmdup -t 4 --prefix= $OUTDIR/$probe$SUF.fa   $OUTDIR/$probe$SUF.fa 
	mv  ${probe}$SUF.rmdup*  "$OUTDIR/"  		#### bo opcja --prefix nie działa
}
rmdup $CONTR $SUF &
rmdup $TREAT $SUF &
wait
DUPSUF=$SUF.rmdup
echo `date` 'rmdup'

########### Merging treated and control, remove duplicates ##########
PREFIX=merged
#####sga merge usuwa zawarte sekwencje. Dla wersji z trimingiem odczytow trzeba je przygotować inaczej (patrz poniżej -->). Dla odczytów równej długości to jest chyba szybsze.
sga merge -p $OUTDIR/$PREFIX$SUF -t 8 $OUTDIR/${CONTR}${DUPSUF}.fa $OUTDIR/${TREAT}${DUPSUF}.fa
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

#sga rmdup -t 8 --prefix=$OUTDIR/$PREFIX$SUF $OUTDIR/$PREFIX$SUF.fa  
echo `date` 'rmdup merged'
SUF=$SUF.rmdup
mv  $PREFIX$SUF*  $OUTDIR/  		#### bo opcja --prefix nie działa


#### finding which duplicate duplicates which sequence using overlap graph build
sga overlap -m $OVERLAP -t 16 --target-file=$OUTDIR/${PREFIX}$SUF.fa $OUTDIR/${PREFIX}$SUF.dups.fa
echo `date` 'overlap duplicates' 
mv ${PREFIX}$SUF* $OUTDIR/				#### bo opcja --prefix nie działa
DUPSGRAPH=$OUTDIR/${PREFIX}$SUF.dups.${PREFIX}$SUF.asqg
gunzip -c ${DUPSGRAPH}.gz > $DUPSGRAPH



########## Overlaps ###########

#sga overlap -m $OVERLAP -t 16 --prefix=$OUTDIR/${PREFIX}$SUF $OUTDIR/${PREFIX}$SUF.fa 
echo `date` 'overlap'
mv  ${PREFIX}$SUF.asqg.gz  $OUTDIR/${PREFIX}${SUF}_${OVERLAP}.asqg.gz   		#### bo opcja --prefix nie działa
SUF=${SUF}_${OVERLAP}
#gunzip -c $OUTDIR/${PREFIX}$SUF.asqg.gz > $OUTDIR/${PREFIX}$SUF.asqg  



######### Asemble -- get contigs and create simplified chunk graph ##########
#sga assemble -o $OUTDIR/${PREFIX}$SUF $OUTDIR/${PREFIX}$SUF.asqg.gz
#echo `date` 'assemble'
### rozpakowanie "oczyszczonego" grafu nałożeń - assembly wykonuje kilka rzeczy żeby graf uprościć 
#gunzip -c $OUTDIR/${PREFIX}$SUF-graph.asqg.gz > $OUTDIR/${PREFIX}$SUF-graph.asqg





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



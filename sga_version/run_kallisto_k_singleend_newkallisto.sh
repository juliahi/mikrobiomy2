K=$1 #max 31
GENOMES_FILE=$2
DIR=$3
mkdir -p $DIR


echo $K $GENOMES_FILE $DIR

INDEX_FILE=$DIR/kallisto_index_$K.idx

LOGFILE=$DIR/kallisto_singleend.log

#KALLISTO="/home/julia/kallisto_linux-v0.43.0/kallisto"
KALLISTO="/home/julia/kallisto_kod/src/kallisto"
KALLISTO="/home/julia/kallisto_kod_merge/src/kallisto"



##build index for genome
if [ ! -f "$INDEX_FILE" ]; then
    echo "build index for $K, $INDEX_FILE, $GENOMES_FILE" >> $LOGFILE
    $KALLISTO index -k $K -i $INDEX_FILE $GENOMES_FILE  2>> $LOGFILE
    echo 'index prepared'
fi

DIR=$DIR/singleend_newkallisto
mkdir -p $DIR


#quantify
licznik=0

for file in $INDIR/*depl_1.fq.gz; do
    FILENAME=${file%_1.fq.gz}
    OUTNAME=$DIR/`basename ${FILENAME}`_kallisto_${K}_out
    echo $file
    
    SECOND=$INDIR/`basename ${FILENAME}`_2.fq.gz

    LOGFILE=${OUTNAME}.log
    SAMFILE=${OUTNAME}.sam
    COMMAND1="$KALLISTO quant -i $INDEX_FILE -o ${OUTNAME} --pseudobam --single -l 100 -s 0.2 ${file} $SECOND "
    #$COMMAND1 2>>${LOGFILE} &
    licznik=$((licznik+1))
    if [ $licznik -eq 5 ]; then
        wait
        licznik=0
    fi
done
wait


licznik=0
for file in $INDIR/*depl_1.fq.gz; do
	FILENAME=${file%_1.fq.gz}
	
	OUTNAME=$DIR/`basename ${FILENAME}`_kallisto_${K}_out
    	(  samtools sort -@ 2 ${OUTNAME}/pseudoalignments.bam -o ${OUTNAME}_pseudoal_sorted.bam && samtools index ${OUTNAME}_pseudoal_sorted.bam ) &
		
	licznik=$((licznik + 1))
    	if [ $licznik -eq 5 ]; then
   		wait
        	licznik=0
    	fi
done
wait



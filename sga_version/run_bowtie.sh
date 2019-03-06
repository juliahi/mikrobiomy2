#WYMAGANIA
# biopython
# ht-seq

#INDIR="/mnt/chr4/mikrobiomy-2/Wyniki_sekwencjonowania/demultiplexed"
INDIR="/home/julia/Wyniki_sekwencjonowania"
export INDIR


GENOMES_FILE=$1
OUTDIR=$2/bowtie_results
mkdir -p $OUTDIR


### map with bowtie2
BOWTIEDIR=/home/julia/lib/bowtie2-2.2.9

INDEXFILE=$OUTDIR/bowtie2_index
LOGFILE=$OUTDIR/bowtie2_sensitive.log
$BOWTIEDIR/bowtie2-build --threads 8 $GENOMES_FILE $INDEXFILE 2>> $LOGFILE >> $LOGFILE
echo "index build"


for file in $INDIR/*depl_*1.fq.gz; do
    FILENAME=${file%_1.fq.gz}
    probe=`basename $FILENAME`
    OUTNAME=$OUTDIR/${probe}_bowtie

    echo $probe 
    SAMFILE=${OUTNAME}.sam
    $BOWTIEDIR/bowtie2 --threads 8 --very-sensitive -N 1 -x $INDEXFILE -1 ${FILENAME}_1.fq.gz -2 ${FILENAME}_2.fq.gz -S $SAMFILE 2>> $LOGFILE

    samtools view -b ${SAMFILE} > ${OUTNAME}.bam 
    rm ${SAMFILE}
    samtools sort -@ 2 ${OUTNAME}.bam -o ${OUTNAME}_sorted.bam
    samtools index ${OUTNAME}_sorted.bam
done

wait
python parse_bowtie_result.py $LOGFILE






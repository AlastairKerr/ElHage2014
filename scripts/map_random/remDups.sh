
java -jar /usr/local/picard-tools/picard-tools-1.107/MarkDuplicates.jar I=$1.bam O=$1.NR.bam AS=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$1.dupmetrics

samtools index $1.NR.bam
samtools flagstat $1.NR.bam > $1.NR.bam.flagstat 

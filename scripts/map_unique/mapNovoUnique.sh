#usage sh mapNovoUnique.sh <input.fastq> <output.prefix>

novoalign -d /homes/genomes/s.cerevisiae/sacCer3/novoindexes/sacCer3.novoindex -f $1 -F STDFQ -o SAM -r None -c 8 > $2.novoUnique.sam

samtools view -@ 10 -Shb $2.novoUnique.sam -F 4 -o $2.novoUnique.bam

samtools sort -@ 10 $2.novoUnique.bam $2.sorted.novoUnique

rm $2.novoUnique.bam $2.novoUnique.sam                                                                                                                     

java -jar /usr/local/picard-tools/picard-tools-1.107/MarkDuplicates.jar I=$2.sorted.novoUnique.bam O=$2.sorted.novoUnique.dups.bam AS=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$2.sorted.novoUnique.dupmetrics

mv $2.sorted.novoUnique.dups.bam $2.sorted.novoUnique.bam

samtools index $2.sorted.novoUnique.bam

samtools flagstat $2.sorted.novoUnique.bam > $2.sorted.novoUnique.bam.flagstat 

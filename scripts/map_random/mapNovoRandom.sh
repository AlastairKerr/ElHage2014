#usage sh mapNovoRandom.sh <input.fastq> <output.prefix>

novoalign -d /homes/genomes/s.cerevisiae/sacCer3/novoindexes/sacCer3.novoindex -f $1 -F STDFQ -o SAM -r Random -c 8 > $2.novoRandom.sam

samtools view -@ 10 -Shb $2.novoRandom.sam -F 4 -o $2.novoRandom.bam

samtools sort -@ 10 $2.novoRandom.bam $2.sorted.novoRandom

rm $2.novoRandom.bam $2.novoRandom.sam                                                                                                                     

java -jar /usr/local/picard-tools/picard-tools-1.107/MarkDuplicates.jar I=$2.sorted.novoRandom.bam O=$2.sorted.novoRandom.dups.bam AS=true VALIDATION_STRINGENCY=LENIENT METRICS_FILE=$2.sorted.novoRandom.dupmetrics

mv $2.sorted.novoRandom.dups.bam $2.sorted.novoRandom.bam

samtools index $2.sorted.novoRandom.bam

samtools flagstat $2.sorted.novoRandom.bam > $2.sorted.novoRandom.bam.flagstat 

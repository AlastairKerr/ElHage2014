

intersectBed -abam $1.bam -b $2 -sorted -v  > $1.blacklist.bam

samtools sort -@ 10 $1.blacklist.bam $1.blacklist.sorted

mv $1.blacklist.sorted.bam $1.blacklist.bam

samtools index $1.blacklist.bam

samtools flagstat $1.blacklist.bam > $1.blacklist.bam.flagstat 

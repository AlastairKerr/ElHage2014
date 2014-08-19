#usage sh raw.sh <directory_of_indexed_bam_files>


for i in `ls $1 | grep -v .bai`
do 


#Deep Tools bamCoverage to convert to bigWig
genomeCoverageBed -ibam $1/$i -bg > $i.bg
#Change chrMito to chrM
sed -i s/chrMito/chrM/ $i.bg
#replace converted bigwig
wigToBigWig $i.bg sacCer3adj.len $i.bw;

rm $i.bg


done

 

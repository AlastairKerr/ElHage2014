
#usage sh intergenicScale.sh <directory_of_indexed_bam_files>

for i in `ls $1 | grep -v .bai`
do 


#Deep Tools bamCoverage to convert to bigWig
bamCoverage -bs 50 -p 15 --smoothLength 100 -b $1/$i --scaleFactor 1 -f 200 --outFileName $i.bw
#Convert to wiggle
bigWigToWig $i.bw $i.wig
#Change chrMito to chrM
sed -i s/chrMito/chrM/ $i.wig
#replace converted bigwig
wigToBigWig $i.wig sacCer3adj.len $i.bw;

#Window genome
perl bigWigAverage.pl --in $i.bw --bed sacCer3.intergenic.bed --out $i.intergenic.mean
#Scale windows to genome wide median
sh make_scaled_wigs.sh $i.wig `cat $i.intergenic.mean`
#remove intermediate files
rm $i.wig

done

 

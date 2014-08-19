#usage sh pipeline.sh <sample_file_of_input_trimmed_fastqs> <output_dir>


for i in `cat $1`; 

do 

#map reads
sh mapNovoUnique.sh $2/$i.trimmed.fq $2/$i
#remove rDNA reads
sh blacklist.sh $2/$i.sorted.novoUnique rDNA.bed
#remove duplicates
sh remDups.sh $2/$i.sorted.novoUnique.blacklist 

done;

 

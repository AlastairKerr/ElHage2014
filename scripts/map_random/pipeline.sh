
#usage sh pipeline.sh <sample_file_of_input_fastq_files> <output_dir>

for i in `cat $1`; 

do 
#trim reads
sh trimmomatic.sh $2/$i.fq $2/$i $2/$i.trim.log; 
#map reads
sh mapNovoRandom.sh $2/$i.trimmed.fq $2/$i
#remove rDNA reads
sh blacklist.sh $2/$i.sorted.novoRandom rDNA.bed
#remove duplicates
sh remDups.sh $2/$i.sorted.novoRandom.blacklist 

done;

 

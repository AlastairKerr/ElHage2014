
java -classpath /opt/software/Trimmomatic/trimmomatic-0.22.jar org.usadellab.trimmomatic.TrimmomaticSE -phred33 -threads 14 -trimlog $3 $1 $2.trimmed.fq ILLUMINACLIP:/homes/swebb/tool_files/adapters_illumina_genomic_w_revcomps.fasta:2:40:15 HEADCROP:1 SLIDINGWINDOW:5:20 MINLEN:30

fastqc $2.trimmed.fq 

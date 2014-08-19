
*fastq files were aligned to sacCer3 genome using mapNovoUnique.sh for unique alignments only and mapNovoRandom.sh to report a single random alignment.

*bigWig files were created and normalised to the intergenic mean using intergenicScale.sh

*bigWigPlotter.pl was used to bin regions of interest in order to create average profile plots

*mergeWindows.pl was used to combine exon1, intron1, exon2 binned regions in to one region

*addWindowFlag.pl and mapColumn.pl were used to add columns to the window files to allow regions of interest to be filtered (i.e. gene type, expression level etc). 

*Average profile plots were created with custom R scripts based on plotWindows.pl

*Annotation files are provided

*Thermo stability plots were created and plotted using ThermoWin.pl

####Dependencies####
Samtools
BedTools
Picard Tools
Novoalign
Fastqc
Trimmomatic
UCSC tools (wigToBigWig, genomeCoverageBed, bigWigAverageOverBed)
Deep Tools
Perl
R


 

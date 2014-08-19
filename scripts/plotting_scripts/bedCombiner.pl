#!/usr/bin/perl -w 

=head1 NAME

mergeIntervals.pl

=head1 SYNOPSIS

Combines interval files by looking for overlapping regions

=head1 Options 

Flagged arguments

=head2 --add (required)

comma separated list of input files

=head2 --outfile (required)

output file

=head2 --chrcol                                                                                                               

Set chromosome column (default=1)
 
=head2 --startcol

Set start column (default=1)                                                                                                                
=head2 --endcol             

Set end column (default=2)
                                                                                                   
=cut

#Combines bed files by looking for overlapping regions

use Getopt::Long;
use strict;

my %opt = (
	chrcol => 1,
	startcol => 2,
	endcol => 3,
	);

my @add;
            
GetOptions(\%opt, "add=s{,}" => \@add, "outfile=s", "chrcol=i","startcol=i","endcol=i");


die "ERROR: must enter at least 1 input file"
    unless (@add > 0);

$opt{chrcol} --;
$opt{startcol} --;
$opt{endcol} --;

my %chromosomes; #represents the chromosomes covered in input files
my %region; #represents each base pair covered in a specific region e.g. 1 chromosome



#
# Read all files and get chromosome numbers
#
my $aLength = @add;
for (my $i=0;$i<$aLength;$i++){
	&getChroms($add[$i]);
}




#
# For each chromosome read all files, do sums and print results
#
open(OUT, ">$opt{outfile}" ) or die "cannot write to $opt{outfile}";	
foreach my $chromo (sort {$a cmp $b} keys %chromosomes){
	for (my $i=0;$i<$aLength;$i++){
		my $addFile = $add[$i];		
		&readBed($addFile, $chromo);
	}

	#print combined chromosome values to output
	&printFile($chromo);

	#reset variables for new chromosome
	%region = ();
}

close OUT;




#read in each line of the file, get the chr numbers represented and add them to the chromosome hash as keys
sub getChroms
{

#    my $file = $_[0];
#    my $col = $opt{chrcol} + 1;
#    system "cut -f $col -d \"	\" $file | uniq > chroms.tmp";

#    open(CHR, "chroms.tmp") or die "cannot create chromosome list";
#    while(<CHR>){
#	chomp;
#	if (!exists $chromosomes{$_}){
#		$chromosomes{$_} = 0;
#	}
#    }
#    close CHR;
#    system "rm chroms.tmp";

    my $file = $_[0];
    open(BED, $file) or die "cannot open $file\n";
    while(<BED>){
	chomp;
	my $c  = (split /\t/)[$opt{chrcol}]; #get chr numbers
	if (!exists $chromosomes{$c}){
		$chromosomes{$c} = 0;
	}
    }
    close BED;
}




# read in each line of the file and update value in a hash with bp co-ords as key depending on combining operation
sub readBed 
{
    my $readFile = $_[0];
    my $chr = $_[1];	
    open(BEDS, $readFile) or die "cannot open $readFile\n";
    while(<BEDS>){
	next if /#/; #skip comment lines	
	chomp;
	my($c, $s, $e) = (split /\t/)[
				      $opt{chrcol},	
				      $opt{startcol}, 
				      $opt{endcol}, 
				      ]; #Chr Start Stop 

	if ($c eq $chr){

		#add bases to region hash
		$region{$s}||=0;
		$region{$s}+=1; #Add a start (1) to the region hash
		$region{$e}||=0;
		$region{$e}-=1; #Add an end (-1) to the region hash

	}
    }
    close BEDS;
}


#find overlapping regions and print in bed format
sub printFile
{
	my $chrom = $_[0]; #chromosome number
	my $contig = 0;	#boolean, 1 = in middle of contig
	my $wigStart = 0; #start of contig
	my $wigEnd=0; #end of contig

	my $contigCount=0; #sums the values at each start (+1) and end (-1). End of contig when contigCount = 0.
		
	#for each start and end value	
	foreach my $base (sort {$a <=> $b} keys %region){	
	    $contigCount += $region{$base}; #sum the values of starts and ends
	    if ($contig == 0){  #if 0 then new contig
		$wigStart = $base;  #set start
		$contig = 1;
	    }
		
	    if ($contigCount == 0){ #end of contig
		$wigEnd = $base;  #set end
		print OUT ("$chrom\t$wigStart\t$wigEnd\n");
		#reset for next contig
	        $contig = 0;
	    }	
	       
	}
}


 

#!/usr/bin/perl -w



my    %dna_dna = (	'AA' => [7.9, 0.0222],
	'TT' => [7.9, 0.0222],
	'AT' => [7.2, 0.0204],
	'TA' => [7.2, 0.0213],
	'CA' => [8.5, 0.0227],
	'TG' => [8.5, 0.0227],
	'GT' => [8.4, 0.0224],
	'AC' => [8.4, 0.0224],
	'CT' => [7.8, 0.021],
	'AG' => [7.8, 0.021],
	'GA' => [8.2, 0.0222],
	'TC' => [8.2, 0.0222],
	'CG' => [10.6, 0.0272],
	'GC' => [9.8, 0.0244],
	'GG' => [8.0, 0.0199],
	'CC' => [8.0, 0.0199]
	);
    
## RNA- DNA 
my   %rna_dna = (
	'AA' => [7.8, 0.0219],
	'TT' => [11.5, 0.0364],
	'AT' => [8.3, 0.0239],
	'TA' => [7.8, 0.0232],
	'CA' => [9.0, 0.0261],
	'TG' => [10.4, 0.0284],
	'GT' => [7.8, 0.0216],
	'AC' => [5.9, 0.0123],
	'CT' => [7.0, 0.0197],
	'AG' => [9.1, 0.0235],
	'GA' => [5.5, 0.0135],
	'TC' => [8.6, 0.0229],
	'CG' => [16.3, 0.0471],
	'GC' => [8.0, 0.0171],
	'GG' => [12.8, 0.0319],
	'CC' => [9.3, 0.0232]
	);
use Getopt::Long;
use Bio::Seq;
use Bio::SeqIO;
use strict;

my %opt = (
	outfile => "Thermo", 
	window  =>  9,
	step    =>  9,
	);
            



GetOptions(\%opt, 
	   "infile=s", 
	   "outfile=s", 
	   "window=i", 
	   "step=i", 
	   );


die "ERROR: must use '--infile' to set name of your sequence file"
  if !$opt{infile};

my $in  = Bio::SeqIO->new('-file' => $opt{infile}, 
                        '-format' => 'Fasta'); 

open (OUTGC, ">$opt{outfile}_gc.bg")|| die "cannot write to  $opt{outfile}";
open (OUTD, ">$opt{outfile}_dnadna.bg")|| die "cannot write to  $opt{outfile}";
open (OUTR1, ">$opt{outfile}+rnddna.bg")|| die "cannot write to  $opt{outfile}";
open (OUTR2, ">$opt{outfile}-rnddna.bg")|| die "cannot write to  $opt{outfile}";
open (OUTR1R, ">$opt{outfile}+rnddna_ratio.bg")|| die "cannot write to  $opt{outfile}";
open (OUTR2R, ">$opt{outfile}-rnddna_ratio.bg")|| die "cannot write to  $opt{outfile}";


while (my $seqobj = $in->next_seq()){
    my $contig_len = $seqobj->length;
    my $start = 1;
    while($start + $opt{window} < $contig_len){
	my $seq = $seqobj->subseq($start, $start+ $opt{window}-1);
	next if (length($seq) < 1);
	
	my $n = $seq =~ tr/N/N/; #count Ns

	if ($n > length($seq)/3){
	    $start += $opt{step};
	    next;
	}

#get reverse complement of window
	my $revseq = reverse($seq);
	$revseq =~ tr/ATGC/TACG/;


#	my $o = observed_CG($seq);
	my($e, $g_c) = expected_CG($seq);
	my $dd_G = thermo_dnadna($seq);
	my $rd_G1 = thermo_rnadna($seq);
	my $rd_G2 = thermo_rnadna($revseq);

	print OUTGC join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 $g_c
	    ),"\n";

	print OUTD join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 $dd_G
	    ),"\n";

	print OUTR1 join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 $rd_G1
	    ),"\n";

	print OUTR2 join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 $rd_G2
	    ),"\n";

	print OUTR1R join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 sprintf("%.3f", $rd_G1 / $dd_G)
	    ),"\n";
        

	print OUTR2R join("\t", 
			 $seqobj->display_id,
			 $start,
			 $start + $opt{window},
			 sprintf("%.3f", $rd_G2 / $dd_G)
	    ),"\n";

	$start += $opt{step};
    }
}


#sub observed_CG{
#    my $seq = shift;
#    my $CG = 0;
#    $CG += $seq =~ s/CG/CG/g;
#    $CG;
#}



sub expected_CG{
    my $CG = 0;
    my $seq = shift;
    my $G = $seq =~ tr/G/G/;
    my $C = $seq =~ tr/C/C/;
    my $len = length($seq);
    return("error") unless ($len);
    $CG = ( ($G/$len) * ($C/$len) * ($len-1));

    return($CG,  sprintf("%.3f", (($G + $C)/$len) ));
}


    


sub thermo_rnadna{
    my $sub_seq = shift;
    my $H = 0;
    my $S = 0;
    my $G = 0;
    for (my $j = 0; $j < length($sub_seq)-1; $j++)
    {
	my $pair = substr($sub_seq,$j,2);
	$pair = uc($pair);              # convert to upper case
        unless($pair =~ /nx- /i){        
	    $H += $rna_dna{$pair}->[0];
	    $S += $rna_dna{$pair}->[1]; 
	}
    }
    if ( $S != 0 ) { # delta G 
	$G =sprintf("%.3f", $H*(1-310/($H/$S)));
    }else{
	$G="NA";;
    }
    $G;
}

  
sub thermo_dnadna{
    my $sub_seq = shift;
    my $H = 0;
    my $S = 0;
    my $G = 0;
    for (my $j = 0; $j < length($sub_seq)-1; $j++)
    {
	my $pair = substr($sub_seq,$j,2);
	$pair = uc($pair);              # convert to upper case
        unless($pair =~ /nx- /i){        
	    $H += $dna_dna{$pair}->[0];
	    $S += $dna_dna{$pair}->[1]; 
	}
    }
    if ( $S != 0 ) { # delta G 
	$G =sprintf("%.3f", $H*(1-310/($H/$S)));
    }else{
	$G="NA";;
    }
    $G;
}

  

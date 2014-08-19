#!/usr/bin/perl -w

=head1 NAME

bedExons.pl

=head1 SYNOPSIS

Splits Bed 12 file in to tabular file with position of exon1,intron1,exon2....

=head1 Options

Flagged arguments

=head2 --in (required)

name of the input file

=head2 --out (required)

name of the output file

=cut

use strict;
use Getopt::Long;

my %opt = ();
            
GetOptions(\%opt, "in=s", "out=s");

open (IN, "$opt{in}") || die "cannot open $opt{in}";
open (OUT, ">$opt{out}") || die "Could not write to $opt{out}";


while (<IN>){
    chomp;
    next if /^#/;
    my @bed = split(/\t/,$_);
    my $l = @bed;
    if ($l < 12){
	print STDERR "Please use 12 column bed file with exon information";
	exit;
    }
    my $strand = $bed[5];
    my @exstart = split(/,/,$bed[11]);
    my @exlen = split(/,/,$bed[10]);
    my $exnum = $bed[9];

    my %ex=();
    my %in=();

    for (my $i=0;$i<$exnum;$i++){
	my $estart = $bed[1] + $exstart[$i];
	my $eend = $bed[1] + $exstart[$i] + $exlen[$i];
        if($strand eq '+'){
	   push (@{$ex{$i+1}},$estart);
	   push (@{$ex{$i+1}},$eend);
	}
	else{
	    push (@{$ex{$exnum-$i}},$estart);
	    push (@{$ex{$exnum-$i}},$eend);
	}
    }
    

    print OUT "$bed[0]\t$bed[1]\t$bed[2]\t$bed[3]\t$bed[4]\t$bed[5]";
    for(my $i=1;$i<=$exnum;$i++){
	print OUT "\t$ex{$i}[0]\t$ex{$i}[1]";
	if ($i < $exnum){
	    if ($strand eq '+'){
		print OUT "\t$ex{$i}[1]\t$ex{$i+1}[0]";
	    }
	    else{
		print OUT "\t$ex{$i+1}[1]\t$ex{$i}[0]";
	    }
	}
    }

    print OUT "\n";

}

close OUT;
 

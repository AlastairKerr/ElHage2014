#!/usr/bin/perl -w

=head1 NAME

bedSpliceSite.pl

=head1 SYNOPSIS

Splits Bed 12 file in to tabular file with position of splice site.

=head1 Options

Flagged arguments

=head2 --in (required)

name of the input file

=head2 --out (required)

name of the output file

=head2 -s (required)

type of site (tss,3,5,tes)

=head2 -n

number of site (default=1)

=cut

use strict;
use Getopt::Long;

my %opt = (n=>1);
            
GetOptions(\%opt, "in=s", "out=s","s=s","n=i");

open (IN, "$opt{in}") || die "cannot open $opt{in}";
open (OUT, ">$opt{out}") || die "Could not write to $opt{out}";


while (<IN>){
    chomp;
    next if /^#/;
    my @bed = split(/\t/,$_);
    my $strand = $bed[5];
    my @exstart = split(/,/,$bed[11]);
    my @exlen = split(/,/,$bed[10]);
    my $exnum = $bed[9];

    my $site;

    if($opt{s} eq "tss"){
	if($strand eq "-"){
	    $site = $bed[2]-1;
	}
	else{
            $site = $bed[1];
	}
    }
    elsif($opt{s} eq "tes"){
	if($strand eq "-"){
            $site = $bed[1];
        }
	else{
            $site = $bed[2]-1;
        }
    }
    elsif($opt{s} eq "3"){
	if($opt{n}>=$exnum){next;}
	if($strand eq "-"){
            $site = $bed[1]+$exstart[$exnum-$opt{n}-1]+$exlen[$exnum-$opt{n}-1]-1;
        }
        else{
            $site = $bed[1]+$exstart[$opt{n}];
        }
    }
    elsif($opt{s} eq "5"){
	if($opt{n}>=$exnum){next;}
        if($strand eq "+"){
            $site = $bed[1]+$exstart[$opt{n}-1]+$exlen[$opt{n}-1]-1;
        }
        else{
            $site = $bed[1]+$exstart[$exnum-$opt{n}];
        }
    }


    if($strand eq "+"){
	my $end = $site +1;
	print OUT "$bed[0]\t$site\t$end\t$bed[3]\t$bed[4]\t$bed[5]\n";
    }
    else{
	my $end = $site +1;
        print OUT "$bed[0]\t$site\t$end\t$bed[3]\t$bed[4]\t$bed[5]\n";
    }
}

close OUT;
 

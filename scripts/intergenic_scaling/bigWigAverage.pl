#!/usr/bin/perl -w 


use Getopt::Long;
use strict;

#get user options
my %opt = ( );
            
GetOptions(\%opt, "in=s", "bed=s", "out=s");


open(OUT, ">$opt{out}");

system "bigWigAverageOverBed $opt{in} $opt{bed} $opt{in}.intergenic.tab";

open(MEAN, "$opt{in}.intergenic.tab");

my $s = 0;
my $c = 0;

while(<MEAN>){
    chomp();
    my $cov = (split /\t/,$_)[2];
    my $sum = (split /\t/,$_)[3];

    $s += $sum;
    $c += $cov;
}

my $mean = $s/$c;

print OUT "$mean";





 

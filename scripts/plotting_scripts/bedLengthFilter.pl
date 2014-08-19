#!/usr/bin/perl -w
=head1 NAME

bedLengthFilter.pl

=head1 SYNOPSIS

Splits a bed file in two based on the length of the intervals: ge = greater or equal to, lt = less than.

=head1 Options 

Flagged arguments

=head2 --in (required)

input file

=head2 --lt (default=lt.100.bed)

lt output file

=head2 --ge (default=ge.100.bed)

ge output file

=head2 -l

length of region to filter on (default=100) 

=cut
use strict;
use Getopt::Long;

my %opt = (l=>100,);
            
GetOptions(\%opt, "in=s","lt=s", "l=i", "ge=s");

open(IN, $opt{in}) || die "Cannot open $opt{in}";

if(!exists $opt{lt}){
    $opt{lt} ="lt.$opt{l}.bed";
}
if(!exists $opt{ge}){
    $opt{ge} ="ge.$opt{l}.bed";
}

open(SHORT,">$opt{lt}") || die "Cannot write to $opt{lt}";
open(LONG,">$opt{ge}") || die "Cannot write to $opt{ge}";


while (<IN>){
    chomp($_);
    my ($s,$e) = (split /\t/,$_)[1,2];
    my $length = $e-$s;

    if($length<$opt{l}){
	print SHORT "$_\n";
    }
    else{
	print LONG "$_\n";
    }
}

close IN;
close SHORT;
close LONG;
 

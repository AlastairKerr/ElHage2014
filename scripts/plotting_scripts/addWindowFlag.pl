#!/usr/bin/perl -w

=head1 NAME

addWindowFlag.pl

=head1 SYNOPSIS

Add an extra column to flag certain regions for downstream filtering (e.g. RP genes).

=head1 Options 

Flagged arguments

=head2 --in (required)

Input window file

=head2 --out (required)

name of the output file

=head2 -f

Character to add as a flag

=head2 -af

Character to add as anti-flag (default=NA)

=head2 --id

Column 1 = ID to flag, column 2 = flag (if -f not set)

=head2 --bed

Use a bed file to select flagged IDs

=cut

use strict;
use Getopt::Long;

my %opt = ();
GetOptions(\%opt, "in=s", "out=s","f=s","af=s","bed=s","id=s");

open(OUT,">$opt{out}") || die "cannot write to $opt{out}";
open(IN,$opt{in}) || die "cannot open $opt{in}";

#Get IDs to flag
my %ref=();
if(exists $opt{bed}){
    open(BED,$opt{bed}) || die "cannot open $opt{bed}";
    while(<BED>){
	chomp();
	my $id = (split /\t/,$_)[3];
	if (!defined $id){
	    print STDERR "Bed file not in correct format";
	    exit;
	}
	$ref{$id}=1;
    }
    close BED;
}
elsif(exists $opt{id}){
    open(ID,$opt{id}) || die "cannot open $opt{id}";
    while(<ID>){
        chomp();
        my ($id,$fl) = (split /\t/,$_)[0,1];
	$ref{$id}=$fl;
    }
    close ID;
}
else{
    print STDERR "Must use --bed or --id options";
    exit;
}


while(<IN>){      
    chomp();
    my @line = (split /\t/,$_);
    my $flag = "";
    my $idc = 3;

    if(exists $ref{$line[$idc]}){
	if(exists $opt{f}){
	    $flag = $opt{f};
	}
	else{
	    $flag = $ref{$line[$idc]};
	}
    }
    elsif(exists $opt{af}){
	$flag = $opt{af};
    }
    else{
	$flag = "NA";
    }
    my $f = 7+($line[6]*3);
    $line[$f] ++;
   
    my $l = @line;
    for(my $i=0;$i<=$f;$i++){
	print OUT "$line[$i]\t";
    }
    print OUT "$flag";
    for(my $i=$f+1;$i<$l;$i++){
	print OUT "\t$line[$i]";
    }
    print OUT "\n";
}
close IN;
close OUT;


 

#!/usr/bin/perl -w

=head1 NAME

mergeWindows.pl

=head1 SYNOPSIS

Merge sliding window files based on bed ID column, each section of windows must belong to a unique ID. Add flags to windows after merging or they will be lost.

=head1 Options 

Flagged arguments

=head2 --in (required)

Comma separated list of files in order of merge

=head2 --out (required)

name of the output file

=cut

use strict;
use Getopt::Long;

my %opt = ();
my @in = ();            
GetOptions(\%opt, "in=s{,}" => \@in, "out=s");

#number of files
my $n = @in;

my %wins=();
my %reg=();
my %count=();
my %bed=();


for(my $i=0; $i<$n;$i++){
 
   open(IN,$in[$i]) || die "cannot open $in[$i]";
   while(<IN>){      
      chomp();
      my @line = (split /\t/,$_);

      my $f = 7+($line[6]*3); #flag

      my $c = $f+$line[$f]+1;
      my $id = $line[3];

      my $l = @line - $c;
      my $reg = "$line[1]\t$line[2]\t$l";

      my @wins = @line[$c..$#line];
      my @bed = @line[0..5];


      if($i==0){
	  @{$bed{$id}} = @bed;
	  $reg{$id} = $reg;
	  $wins{$id}= join("\t",@wins);
	  $count{$id} = 1;
      }
      elsif(exists $reg{$id}){
	  $reg{$id}.="\t$reg";
	  $wins{$id}.="\t".join("\t",@wins);
	  $count{$id} ++;
      }
   }
   close(IN);
}

open(OUT,">$opt{out}") || die "cannot write to $opt{out}";

foreach my $id (keys %count){
    if($count{$id}<$n){
	print STDERR "$id does not exist in all files\n";
    }
}

foreach my $id (sort keys %reg){
    my @r = (split /\t/,$reg{$id});
    my ($s,$e);
    my $l = @r;
 
    if($bed{$id}[5] eq '+'){
	$s = $r[0];
	$e = $r[$l-2];
    }
    else{
        $s = $r[$l-3];
        $e = $r[1];
    }
    $bed{$id}[1] = $s;
    $bed{$id}[2] = $e;
    my $bed = join("\t",@{$bed{$id}});
   
    print OUT "$bed\t$n\t$reg{$id}\t0\t$wins{$id}\n";
}
close OUT;
 

#!/usr/bin/perl -w 

##adds a flank of specified length to interval co-ords

use Getopt::Long;
use strict;

#get user options
my %opt = (
	out => "flanks.interval", 
	flanks => 0,
	flanksOnly => 0,
	returntype => 3,
	chrom => 1,
	start => 2,
	end => 3,
	strand => -1,
	);

GetOptions(\%opt, "intervals=s", "flanks=i", "flanksOnly=i", "chrom=i", "start=i", "end=i", "returntype=i", "strand=i", "out=s" );

my $strand;
if ($opt{strand}<=0){
	$strand = "+";
}


open(INTS, "$opt{intervals}" ) or die "cannot open $opt{ints}";
open(FLANKS, ">$opt{out}" ) or die "cannot create flanks file";
	while(<INTS>){
		next if m/#/;
		chomp;
		my @line = (split /\t/);
		my $c = $line[$opt{chrom}-1];
		my $start=$line[$opt{start}-1];
		my $end=$line[$opt{end}-1];
		if ($opt{strand}>= 1){
			$strand = $line[$opt{strand}-1];
		}
		die "Strand notation must be + or -\n" if ($strand ne "+" and $strand ne "-");

######enter strandedness loops here:

		if ($strand eq "+"){

			if ($opt{returntype} == 1 || $opt{returntype} == 3){ #flank at start or both
				$start = $line[$opt{start}-1] - $opt{flanks};
				$start = 0 if ($start<0);
			}
			if ($opt{returntype} == 2 || $opt{returntype} == 3){ #flank at end or both
				$end = $line[$opt{end}-1] + $opt{flanks};
			}
			my $others = "";
			my $name = "";
			my $length = @line;
			for (my $i=0;$i<$length;$i++){
				if ($i!=$opt{chrom}-1 && $i!=$opt{start}-1 && $i!=$opt{end}-1){
					$others = $others."\t$line[$i]";
				}
			}
		
			if ($opt{flanksOnly} == 1){
				if ($opt{returntype} == 1){
					$end = $line[$opt{start}-1];
				}
				elsif ($opt{returntype} == 2){
					$start = $line[$opt{end}-1];
				}
				elsif ($opt{returntype} == 3){
					my $newEnd = $line[$opt{start}-1];
					my $tempName = "up";
					print FLANKS ("$c\t$start\t$newEnd\t$tempName"."$others\n");
					$start = $line[$opt{end}-1];
					$name = "down";
				}
			}
			if ($opt{flanksOnly} == 1 && $opt{returntype} ==3){
				print FLANKS ("$c\t$start\t$end\t$name"."$others\n");
			}
			else{
				print FLANKS ("$c\t$start\t$end"."$others\n");
			}

		}
		elsif ($strand eq "-"){

			if ($opt{returntype} == 1 || $opt{returntype} == 3){ #flank at start or both
				$end = $line[$opt{end}-1] + $opt{flanks};
			}
			if ($opt{returntype} == 2 || $opt{returntype} == 3){ #flank at end or both
				$start = $line[$opt{start}-1] - $opt{flanks};
				$start = 0 if ($start<0);
			}
			my $others = "";
			my $name = "";
			my $length = @line;
			for (my $i=0;$i<$length;$i++){
				if ($i!=$opt{chrom}-1 && $i!=$opt{start}-1 && $i!=$opt{end}-1){
					$others = $others."\t$line[$i]";
				}
			}

			if ($opt{flanksOnly} == 1){
				if ($opt{returntype} == 1){
					$start = $line[$opt{end}-1];
				}
				elsif ($opt{returntype} == 2){
					$end = $line[$opt{start}-1];
				}
				elsif ($opt{returntype} == 3){
					my $newEnd = $line[$opt{start}-1];
					my $tempName = "down";
					print FLANKS ("$c\t$start\t$newEnd\t$tempName"."$others\n");
					$start = $line[$opt{end}-1];
					$name = "up";
				}
			}
			if ($opt{flanksOnly} == 1 && $opt{returntype} ==3){
				print FLANKS ("$c\t$start\t$end\t$name"."$others\n");
			}
			else{
				print FLANKS ("$c\t$start\t$end"."$others\n");
			}

		}
	}

close INTS;
close FLANKS	


 

#!/usr/bin/perl -w 

=head1 NAME

getRegions.pl

=head1 SYNOPSIS

Create genomic intervals based around entries from a bed file

=head1 Options 

Flagged arguments

=head2 --in (required)

input interval file

=head2 --out (required)

output file

=head2 --sbase

Starting base

=head2 --sreg

Start around which feature (start, end, middle)

=head2 --sdir

direction of start from region (up=upstream, down=downstream)

=head2 --ebase

ending base

=head2 --ereg

End around which feature (start, end, mid)

=head2 --edir

direction of end from region (upstream,downstream)

=head2 --chromInfo

File containing lengths of chromosomes

=head2 --remove

remove regions that extend off the ends of chromosome (default is to truncate regions)

=head2 --strand

Column of strand information if interval features are strand dependent

=cut

##adds a flank of specified length to interval co-ords

use Getopt::Long;
use strict;

#get user options
my %opt = (
	chrom => 1,
	start => 2,
	end => 3,
	strand => -1,
	);

GetOptions(\%opt, "in=s", "sbase=i","sreg=s","sdir=s", "ebase=i","ereg=s","edir=s", "chrom=i", "start=i", "end=i", "strand=i", "out=s", "chromInfo=s", "remove");


my $strand;
if ($opt{strand}<=0){
	$strand = "+";
}

my %lens;


open(LENS,"$opt{chromInfo}") or die "cannot open $opt{chromInfo}";
while(<LENS>){	
	chomp;
	my ($chr,$len) = (split /\t/)[0,1];
	$lens{$chr} = $len;
}



open(INTS, "$opt{in}" ) or die "cannot open $opt{in}";
open(REGS, ">$opt{out}" ) or die "cannot create flanks file";
my $count = 0;
	while(<INTS>){
		$count ++;
		next if m/#/;
		chomp;
		my @line = (split /\t/);
		my $c = $line[$opt{chrom}-1];
		my $start=$line[$opt{start}-1];
		my $end=$line[$opt{end}-1];
		if ($opt{strand}>= 1){
			$strand = $line[$opt{strand}-1];
		}
		die "Start co-ordinate is greater than end co-ordinate at line $count" if ($start >= $end);
		die "Strand notation must be + or -\n" if ($strand ne "+" and $strand ne "-");
		my $mid = int($start + (($end -$start)/2));

######enter strandedness loops here:

		my ($s,$e);
		my $fault = 0;
		if ($strand eq "+"){

			if ($opt{sreg} eq "start"){
				if ($opt{sdir} eq "up"){
					$s =  $start - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$s =  $start + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
					
				}
			}
			elsif ($opt{sreg} eq "end"){
				if ($opt{sdir} eq "up"){
					$s =  $end - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$s =  $end + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
			}
			elsif ($opt{sreg} eq "mid"){
				if ($opt{sdir} eq "up"){
					$s =  $mid - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$s =  $mid + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
			}

			if ($opt{ereg} eq "start"){
				if ($opt{edir} eq "up"){
					$e =  $start - $opt{ebase};
					if ($e<0){
						$e =  0;;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$e =  $start + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
			}
			elsif ($opt{ereg} eq "end"){
				if ($opt{edir} eq "up"){
					$e =  $end - $opt{ebase};
					if ($e<0){
						$e =  0;;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$e =  $end + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
			}
			elsif ($opt{ereg} eq "mid"){
				if ($opt{edir} eq "up"){
					$e =  $mid - $opt{ebase};
					if ($e<0){
						$e =  0;;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
				else{
					$e =  $mid + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
			}

			# if co-ords are the same add 1 to end to represent a single base
			if ($s==$e){
				$e ++;
			}	

			# swap co-ords around if they crossover
			if ($s>$e){
				my $tmp = $s;
				$s = $e;
				$e = $tmp;
			}

			my $others = "";
			my $name = "";
			my $length = @line;
			for (my $i=0;$i<$length;$i++){
				if ($i!=$opt{chrom}-1 && $i!=$opt{start}-1 && $i!=$opt{end}-1){
					$others = $others."\t$line[$i]";
				}
			}
		
			unless (exists $opt{remove} and $fault==1){
				print REGS ("$c\t$s\t$e"."$others\n");
			}


		}


		if ($strand eq "-"){

			if ($opt{sreg} eq "start"){
				if ($opt{sdir} eq "up"){
					$s =  $end + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$s =  $end - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}
			elsif ($opt{sreg} eq "end"){
				if ($opt{sdir} eq "up"){
					$s =  $start + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$s =  $start - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}
			elsif ($opt{sreg} eq "mid"){
				if ($opt{sdir} eq "up"){
					$s =  $mid + $opt{sbase};
					if (exists $lens{$c}){
					if ($s > $lens{$c}){
						$s = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$s =  $mid - $opt{sbase};
					if ($s<0){
						$s =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}

			if ($opt{ereg} eq "start"){
				if ($opt{edir} eq "up"){
					$e =  $end + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$e =  $end - $opt{ebase};
					if ($e<0){
						$e =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}
			elsif ($opt{ereg} eq "end"){
				if ($opt{edir} eq "up"){
					$e =  $start + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$e =  $start - $opt{ebase};
					if ($e<0){
						$e =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}
			elsif ($opt{ereg} eq "mid"){
				if ($opt{edir} eq "up"){
					$e =  $mid + $opt{ebase};
					if (exists $lens{$c}){
					if ($e > $lens{$c}){
						$e = $lens{$c};
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 	
					}
					}
				}
				else{
					$e =  $mid - $opt{ebase};
					if ($e<0){
						$e =  0;
						$fault = 1 ;
						print "Region offset from $c\t$start\t$end\t$strand is beyond end of chromosome\n"; 
					}
				}
			}

			# if co-ords are the same add 1 to end to represent a single base
			if ($s==$e){
				$e ++;
			}	

			# swap co-ords around if they crossover
			if ($s>$e){
				my $tmp = $s;
				$s = $e;
				$e = $tmp;
			}

			my $others = "";
			my $name = "";
			my $length = @line;
			for (my $i=0;$i<$length;$i++){
				if ($i!=$opt{chrom}-1 && $i!=$opt{start}-1 && $i!=$opt{end}-1){
					$others = $others."\t$line[$i]";
				}
			}

			unless (exists $opt{remove} and $fault==1){
				print REGS ("$c\t$s\t$e"."$others\n");
			}

		}
	}

close INTS;
close REGS	


 

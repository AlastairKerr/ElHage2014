#!/usr/bin/perl -w

=head1 NAME

bigWigPlotter.pl

=head1 SYNOPSIS

Performs windowing functions on bigWig file over regions specified in bed file

=head1 Options 

Flagged arguments

=head2 --bw (required)

name of the input bigwig file

=head2 --out (required)

name of the output file

=head2 --bed (required)
                                                                                                                  
name of the bed file containing regions to cover

=head2 -w

number of windows to use (default behaviour = 50 windows)

=head2 --win

Specify a window length instead

=head2 --step

Specify a window step for overlapping windows

=head2 --nodx

Do not reverse windows on - strand regions

=cut

use strict;
use Getopt::Long;
use POSIX qw/ceil/;
use POSIX qw/floor/;


my %opt = (w=>50,strandcol=>6);
            
GetOptions(\%opt, "bed=s", "bw=s","w=i","nodx","win=i","step=i", "out=s","strandcol=i");

open(BED,$opt{bed}) || die "Cannot open $opt{bed}";

my %regions;

#warning that win and step do not reach the end of each region nicely                                                                                   
my $warning = 0;
my $rmdr;
my $first = 0;
my $pid = $$;

#Open output file                                                                                                                                       
open(OUT,">$opt{out}") || die "Cannot write to output file $opt{out}";


#Run with specified window length - Same sized regions only
if(exists $opt{win}){

    #length of the bed region 
    my $len = 0;
    my $count = 1;
    while(<BED>){
	chomp($_);
	my $line = $_;
	my ($c,$s,$e,$n,$z,$st) = (split /\t/, $_)[0,1,2,3,4,$opt{strandcol}-1]; #chrom,start,end and strand
        if(!defined $st){
            $st="+";
        }
        if(!defined $n){
            $n=$count;
            $count ++;
        }
        if(!defined $z){
            $z=0;
        }


        #check bed regions are same size
	my $l = $e - $s;
	if ($len == 0){
          $len = $l;	    
        }
	else{
	    if($l != $len){
		print STDERR "Input regions must be the same length when using parameter --win\n";
		exit;
	    }
	}
	
	#Set the step to window length if not provided
	my $step = $opt{win};
	if(exists $opt{step}){
	    $step = $opt{step};
	}
	if($opt{win} > $l){
	    print STDERR "Window length larger than region size, use a smaller window size\n"; #exit if window is longer than regions
	    exit;
	}

	#calculate number of windows required
	my $num = (int(($l-$opt{win})/$step))+1;
	#check that window and step fit nicely in to reion size
	if($first==0){
	    $rmdr = (($l-$opt{win})%$step);
	    if($rmdr>0){
		$warning=1; #win and step don't fit region size nicely
	    }
	    $first=1;
	}
	
	#Create temp file containing all windows for region as a bed file
       	open(TMP,">$pid.tmp.bed")||die "cannot open temp file";
	my $count=1;
	for(my $i=0;$i<$num;$i++){		
	    my $start = $s + ($i*$step);
	    my $end = $start + $opt{win};

	    print TMP "$c\t$start\t$end\t$count\t0\t$st\n";
	    $count++;
	}

	#Run bigWigAverageOverBed
	system "bigWigAverageOverBed $opt{bw} $pid.tmp.bed $pid.tmp.tab";  

	#Open output and push mean0 column in to an array for each window
	open(BW,"$pid.tmp.tab")||die "cannot open intermediary file $pid.tmp.tab";
        my @array=();
	while(<BW>){
	    chomp($_);
	    my $mean = sprintf("%.2f",(split /\t/, $_)[5]);
	    push (@array,$mean);
	}    
	close BW;

	#Reverse array for regions on minus strand
	if ($st eq "-" and !exists $opt{nodx}){
            @array = reverse @array;
        }

	#print output for that region
	my $windim;
	if ($st eq "-" and !exists $opt{nodx}){
	    $windim = "1\t$e\t$s\t$num\t0";
	}
	else{
	    $windim = "1\t$s\t$e\t$num\t0";
	}
        print OUT "$c\t$s\t$e\t$n\t$z\t$st\t$windim";
	foreach(@array){
            print OUT "\t$_";
	}
        print OUT "\n";


    }	
    #close files
    close BED;
    close OUT;
    if($warning==1){
	print "WARNING: window and step size do not fit nicely in to region size. Final window fell $rmdr bases short of region.\n";
    }
    system "rm $pid.tmp.bed $pid.tmp.tab";
}


elsif(exists $opt{w}){
    my $count =1;
    while(<BED>){
	chomp($_);
	my $line = $_;
	my ($c,$s,$e,$n,$z,$st) = (split /\t/, $_)[0,1,2,3,4,$opt{strandcol}-1];
	if(!defined $st){
            $st="+";
        }
	if(!defined $n){
            $n=$count;
	    $count ++;
        }
        if(!defined $z){
            $z=0;
        }


        #calculate the length of windows required
	my $l = $e-$s;
	my $wl = $l/$opt{w};
	
	#Create temp file containing all windows for region as a bed file 
        open(TMP,">$pid.tmp.bed")||die "cannot open temp file";
        my $count=1;
        for(my $i=0;$i<$opt{w};$i++){
            my $start = $s + ($i*$wl);
            my $end = $start + $wl;
	    my $ws = floor($start);
	    my $we = ceil($end);
	    
            print TMP "$c\t$ws\t$we\t$count\t0\t$st\n";
            $count++;
        }

        #Run bigWigAverageOverBed                                                                                                                          
        system "bigWigAverageOverBed $opt{bw} $pid.tmp.bed $pid.tmp.tab";

	#Open output and push mean0 column in to an array for each window                                                                                  
        open(BW,"$pid.tmp.tab")||die "cannot open intermediary file $opt{bw}.$opt{bed}.tab";
        my @array=();
        while(<BW>){
            chomp($_);
            my $mean = sprintf("%.2f",(split /\t/, $_)[5]);
            push (@array,$mean);
	}
        close BW;

        #Reverse array for regions on minus strand                                                                                                         
        if ($st eq "-" and !exists $opt{nodx}){
            @array = reverse @array;
        }

        #print output for that region
	my $windim;
        if ($st eq "-" and !exists $opt{nodx}){
            $windim = "1\t$e\t$s\t$opt{w}\t0";
	}
	else{
	    $windim = "1\t$s\t$e\t$opt{w}\t0";
	}
        print OUT "$c\t$s\t$e\t$n\t$z\t$st\t$windim";
	foreach(@array){
	    print OUT "\t$_";
	}
	print OUT "\n";

    }
    #close files                                                                                  
    close BED;
    close OUT;
    system "rm $pid.tmp.bed $pid.tmp.tab";
}



 

#!/usr/bin/perl -w

=head1 NAME

mapColumn.pl

=head1 SYNOPSIS

Add columns to a tab file by mapping to a column in another file

=head1 Options 

Flagged arguments

=head2 --input (required)

input file

=head2 --output (required)

output file

=head2 --query (required)

file to map to

=head2 --incol (required)

comma separated list of columns in input file to use as ID

=head2 --column (required)

comma separated list of columns in query file to map to

=head2 --add (required)

comma separated listed of column numbers to add

=head2 --insert

column to insert additional data (default = end of line)

=head2 --newline

Add multiple matches on new line

=head2 --all

return unmapped columns

=head2 --fill

character to fill empty columns

=cut

use strict;
use Getopt::Long;

my %opt = ();

GetOptions(\%opt, "input=s", "output=s", "query=s", "column=s", "incol=s","add=s","newline","all","fill=s","insert=i");


if(exists $opt{insert}){
    $opt{insert} --;
}

open(IN, $opt{input}) or die "cannot open $opt{input}\n";
open(QUERY,"$opt{query}")|| die $!;
open(OUT, ">$opt{output}") or die "cannot write to output";


my @cols = &getOrder($opt{column});
my @incols = &getOrder($opt{incol});
my @addcols = &getOrder($opt{add});
my $valLength = @addcols;
my $length = @cols;
my $testLength = @incols;

die  "You must use the same number of columns in each file" if ($length != $testLength);

my %query;
while (<QUERY>){
    next if /^#/;
    next if /^track/;
    chomp;
    my @line = (split /\t/);
    my $pattern=$line[$cols[0]-1];
    if (!defined $pattern){
	die "Value not present in every column, please use Text Manipulation -> Fill Columns tool"
    }
    for(my $i=1;$i<$length;$i++){
	my $column = (split /\t/)[$cols[$i]-1];
    	if (!defined $column){
		die "Value not present in every column, please use Text Manipulation -> Fill Columns tool"
    	}
	$pattern = "$pattern\t$column"
    }
    my $value=$line[$addcols[0]-1];
    for(my $i=1;$i<$valLength;$i++){
	my $column = (split /\t/)[$addcols[$i]-1];
	$value = "$value\t$column"
    }
    push (@{$query{$pattern}}, $value);
}


while (my $line = <IN>){
    chomp ($line);
    if ($line =~ m/^#/){
	print OUT "$line\n";
	next;
    }
    my $pattern = (split (/\t/,$line))[$incols[0]-1];
    if (!defined $pattern){
	die "Value not present in every column, please use Text Manipulation -> Fill Columns tool"
    }
    for(my $i=1;$i<$length;$i++){
	my $column = (split (/\t/,$line))[$incols[$i]-1];
    	if (!defined $column){
		die "Value not present in every column, please use Text Manipulation -> Fill Columns tool"
    	}
	$pattern = "$pattern\t$column";
    }

    if (exists $query{$pattern}){
	if (exists $opt{newline}){ ##reporting style
		foreach (@{$query{$pattern}}){
		    if(exists $opt{insert}){
			my @line = (split /\t/,$line);
			my $lline = @line;
			if($opt{insert}>$lline){
			    die "Insert column does not exist";
			}
			else{
			    for(my $i=0;$i<$opt{insert};$i++){
				print OUT "$line[$i]\t";
			    }
			    print OUT "$_";
			    for(my $i=$opt{insert};$i<$lline;$i++){
				print OUT "\t$line[$i]";
			    }
			    print OUT "\n";
			}
		    }
		    else{
			print OUT ("$line\t$_\n");
		    }
		}
	}
	else{
	    if(exists $opt{insert}){

		my @line = (split /\t/,$line);
		my $lline = @line;
		if($opt{insert}>$lline){
		    die "Insert column does not exist";
		}
		else{
		    for(my $i=0;$i<$opt{insert};$i++){
			print OUT "$line[$i]\t";
		    }
		    print OUT join("\t",@{$query{$pattern}});
#		    foreach (@{$query{$pattern}}){
#			print OUT ("$_\t");
#		    }
		    for(my $i=$opt{insert};$i<$lline;$i++){
			print OUT "\t$line[$i]";
		    }
		    print OUT "\n";
		}
	    }
	    else{
		print OUT "$line";
		foreach (@{$query{$pattern}}){
			print OUT ("\t$_");
		}
		print OUT "\n";
	    }
	}
    }
    elsif(exists $opt{all}){
	if(exists $opt{insert}){
	    my @line = (split /\t/,$line);
	    my $lline = @line;
	    if($opt{insert}>$lline){
		die "Insert column does not exist";
	    }
	    else{
		for(my $i=0;$i<$opt{insert};$i++){
		    print OUT "$line[$i]";
		}
		for (my $i=0; $i<$valLength; $i++){
		    print OUT "\t";
		    if (exists $opt{fill}){
			print OUT $opt{fill};
		    }
		}
		for(my $i=$opt{insert};$i<$lline;$i++){
		    print OUT "\t$line[$i]";
		}
		print OUT "\n";
	    }
	}
	else{
	    print OUT "$line";
	    for (my $i=0; $i<$valLength; $i++){
		print OUT "\t";	
		if (exists $opt{fill}){
		    print OUT $opt{fill};
		}
	    }
	    print OUT "\n";
	}
    }
}


sub getOrder{
	my $line = shift(@_);
	my @preorder = split (/,/,$line);
	my @order;
	my $count = 0;

	foreach (@preorder){
		if ($_ =~ m/-/){
			my ($start,$end) = (split(/-/,$_))[0,1];
			my $len = $end - $start;
			print STDERR "A negative range $start - $end is specified" if ($len<=0);
			for (my $a = $start; $a<=$end; $a++){
				$order[$count] = $a;
				$count ++;
			}
		}
		else{
			$order[$count] = $_;
			$count ++;
		}
	}
	return @order;
}

 

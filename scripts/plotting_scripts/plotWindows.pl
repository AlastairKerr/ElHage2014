#!/usr/bin/perl -w 

##wrapper for sliding window R plot script

use Getopt::Long;
use strict;

#get user options
my %opt = ();
            
GetOptions(\%opt, "in=s", "box=s","all=s","heat=s","flag=i" );

if(exists $opt{flag}){
    my $c=1;
    my %col=();
    open(F,$opt{in}) || die "Can't open input file";
    open(C,">$opt{in}.tmp.c") || die "Can't write to file ";
    
    while(<F>){
	my $flag=(split /\t/,$_)[$opt{flag}-1];
	if(!exists $col{$flag}){
	    $col{$flag}=$c;
	    $c++;
	}
	print C "$col{$flag}\n";
    }
    close(C);
}


open(R,">$opt{in}.tmp.r") || die "Can't write to file ";


print R "
#R script to create a plots from a table of values
Args<-commandArgs(trailingOnly=T) # get parameters

#print(Args)
# print parameters to screen

table <- read.table(Args[1], header=FALSE, sep=\"\t\",row.names=NULL) #read in the data table
len=dim(table)[2] #get length of table
start=8+(table[1,7]*3)
start=start+table[1,start]+1 #start of data
t<-table[,c(start:len)] #cut out data columns, remove co-ordinate information etc
len2=dim(t)[2]

avg=vector(); #create a list of mean values
for(i in 1:len2){avg<-append(avg,mean(t[,i]))} #for each window, place the mean value in the list

q2=round(len2/4)
q3=round(q2*2)
q4=round(q2*3)


png(Args[2],width = 1000, height = 1000) #open png graphics file
boxplot(t,ylab=\"Score\",xlab=\"Feature length\", varwidth=TRUE,outline=FALSE,col=\"green\",axes=F ) #draw boxplot for each window
axis(2)
axis(1,labels=c(\"0%\",\"25%\",\"50%\",\"75%\",\"100%\"),at=c(0,q2,q3,q4,len2),las=2)
lines(1:len2,avg,col=\"blue\")
dev.off() #close file


library(RColorBrewer)
library(gplots)
hmcol<-brewer.pal(9,\"YlOrRd\")\n";

if(exists $opt{flag}){
print R "colours=read.table(Args[5],header=F)
colours = as.character(colours[,1])";
}

print R "
pdf(Args[3],height=10,width=10)\n";

if(!exists $opt{flag}){
    print R "heatmap.2(as.matrix(t),keysize=1,density.info=\"none\",col=hmcol,trace=\"none\",Colv=NA,symkey=F,labCol=1:len2,labRow=table[,4])\n";
}
else{
    print R "heatmap.2(as.matrix(t),keysize=1,density.info=\"none\",RowSideColors=colours,col=hmcol,trace=\"none\",Colv=NA,symkey=F,labCol=1:len2,labRow=table[,4])\n";
}

print R "dev.off()

png(Args[4],height=1000,width=1000)
g<-dim(t)[1]
max<-max(t)
min<-min(t)
plot(1:len2,t[1,],col=1,ylim=c(min,max),axes=F,type=\"l\",xlab=\"Feature Length\",ylab=\"Score\")
for(i in 2:g){\n";
if(exists $opt{flag}){
    print R " lines(1:len2,t[i,],col=colours[i])\n";
}
else{
    print R "lines(1:len2,t[i,],col=i)\n";
}

print R "}
axis(2)
axis(1,labels=c(\"0%\",\"25%\",\"50%\",\"75%\",\"100%\"),at=c(0,q2,q3,q4,len2),las=2)
dev.off()

";

close(R);

#plot in R                                                                                                                              
if(exists $opt{flag}){
    system "R --slave --file=$opt{in}.tmp.r --args $opt{in} $opt{box} $opt{heat} $opt{all} $opt{in}.tmp.c 2>1";
}
else{
    system "R --slave --file=$opt{in}.tmp.r --args $opt{in} $opt{box} $opt{heat} $opt{all} 2>1";
}


 

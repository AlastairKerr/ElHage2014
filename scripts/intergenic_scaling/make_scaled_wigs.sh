
perl -ne 'next if /#/; chomp; @l = split /\t/;$scaled = sprintf("%.3f",$l[3]/'$2');print "$_\t$scaled\n";' $1 > $1.scaled.tab

cut -f 1,2,3,5 $1.scaled.tab > $1.scaled.wig

NAME=`echo $1 | sed s/.window.tab//`

echo "track type=wiggle_0 name=\"$NAME\" description=\"$NAME\"" > $1.header

cat $1.header $1.scaled.wig > $1.tmp

mv $1.tmp $1.scaled.wig

rm $1.header

wigToBigWig $1.scaled.wig sacCer3adj.len $1.scaled.bw

rm $1.scaled.wig $1.scaled.tab 

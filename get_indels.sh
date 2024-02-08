#!/bin/bash

DY=$1;
WT=$2;
THRESH=2;
CHROM=$3;
cat \
  <(cat  \
    <(tail -n +2 intervals.txt |perl -ane '{print "$F[0]\t",$F[1],"\tstart\n$F[0]\t",$F[2],"\tend\n" if($F[0] eq "'$CHROM'");}')  \
    <(perl -ane '{$h{$F[0]}=1}END{
      open(FILE,"'$DY'");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\t/,$line,2);
        print "$line\n" if(defined($h{$f[0]}));
      }
      }' intervals.txt) \
    |sort -k1,1 -k2,2n -S 10% | \
    awk 'BEGIN{flag=0;}{
      if($3=="start"){
        flag=1;
      }else if($3=="end"){
        flag=0
      };
      if(flag){
        if($1 ~/^#/ || (length($4)!=length($5) && $4!="N" && $5!="N" && $4!~/,/ && $5!~/,/)){
          split($10,a,":");
          if(a[6]>=int("'$THRESH'")) 
            print "DY\t"a[4]"\t"a[6]"\t"$1"\t"$2"\t"$4"\t"$5;
        }
      }
    }' | uniq ) \
  <(cat  \
    <(tail -n +2 intervals.txt |perl -ane '{print "$F[0]\t",$F[1],"\tstart\n$F[0]\t",$F[2],"\tend\n" if($F[0] eq "'$CHROM'");}')  \
    <(perl -ane '{$h{$F[0]}=1}END{
      open(FILE,"'$WT'");
      while($line=<FILE>){
        chomp($line);
        @f=split(/\t/,$line,2);
        print "$line\n" if(defined($h{$f[0]}));
      }
    }' intervals.txt) \
    |sort -k1,1 -k2,2n -S 10% | \
    awk 'BEGIN{flag=0;}{
      if($3=="start"){
        flag=1;
      }else if($3=="end"){
        flag=0;
      };
      if(flag){
        if($1 ~/^#/ || (length($4)!=length($5) && $4!="N" && $5!="N" && $4!~/,/ && $5!~/,/)){
          split($10,a,":");
          if(a[4]>0)
            print "WT\t"a[4]"\t"a[6]"\t"$1"\t"$2"\t"$4"\t"$5;
        }
      }
    }' | uniq ) | \
sort -S 10% -k4,4 -k5,5n -k1,1r| \
uniq -D -f 3 | tee $CHROM.mut.indels.both.txt |\
grep ^DY |\
perl -ane 'BEGIN{$window=750000;$step=50000;}{
  push(@coords,$F[4]);
  push(@freqs,$F[2]/($F[2]+$F[1]));
}
END{
  print "start $coords[0] end $coords[-1]\n";
  for($i=$coords[0]+$window;$i<$coords[-1]-$window;$i+=$step){
    $count=0;
    $value=0;
    for($j=0;$j<=$#coords;$j++){
      if($coords[$j]>=$i-$window && $coords[$j]<$i+$window){
        $count+=1;
        $value+=$freqs[$j];
      }elsif($coords[$j]>=$i+$window){
        if($count >0){
          push(@counts,$value);
        }else{
          push(@counts,$count);
        }
        push(@centers,$i);
        $j=$#coords+1;
      }
    }
  }
  for($i=0;$i<=$#centers;$i++){
    print "$centers[$i] $counts[$i]\n";
  }
}' > $CHROM.mut.indels.index.txt.tmp && \
mv $CHROM.mut.indels.index.txt.tmp  $CHROM.mut.indels.index.txt && \
awk 'BEGIN{max_val=0;max_coord=-1}{if(NF==2){if($2>max_val){max_val=$2;max_coord=$1}}}END{print "Most likely location of the mutation is on chromosome '$CHROM' around position "max_coord}' $CHROM.mut.indels.index.txt && \
awk -F '\t' '{if($1 ~ /^WT/){rowt=$2;aowt=$3}else{if($2<=1 && $3>1) print "DY\t"$2"\t"$3"\tWT\t"rowt"\t"aowt"\t"$4"\t"$5"\t"$6"\t"$7"\t"length($7)-length($6)}}' $CHROM.mut.indels.both.txt >$CHROM.mut.pcr_targets.txt.tmp && \
mv $CHROM.mut.pcr_targets.txt.tmp $CHROM.mut.pcr_targets.txt && rm -f $CHROM.mut.indels.both.txt

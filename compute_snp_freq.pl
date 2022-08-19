#!/usr/bin/env perl
my $ctg="";
my @coords=();
my $window=20000;
my $tol=500;
my $dist=10000;
while($line=<STDIN>){
  chomp($line);
  @f=split(/\s+/,$line);
  process_contig() if(not($f[0] eq $ctg));
  push(@coords,$f[1]);
}
process_contig();

sub process_contig{
    #print "DEBUG: processing contig $ctg with ",$#coords,"\n";
    if($#coords>0){
      #process coords
      my $last_coord=$coords[-1];
      for($j=$window;$j<$last_coord-$window;$j+=$dist){
      my $score=0;
      #print "DEBUG: computing score for $j $last_coord\n";
        for($i=0;$i<$#coords;$i++){
          if($coords[$i]<=$j && $coords[$i+1]>$j){
          #print "DEBUG found center $j index $i value $coords[$i] ",$coords[$i+1],"\n";
          #going backward
            for($k=$i;$k>-1 && $j-$coords[$k]<=$window;$k--){
              #print "DEBUG back $k $score\n";
              $score+=exp(-abs($coords[$k]-$j)/$tol);
            }
          #going forward
            for($k=$i+1;$k<=$#coords && $coords[$k]-$j<=$window;$k++){
              #print "DEBUG forw $k $score\n";
              $score+=exp(-abs($coords[$k]-$j)/$tol);
            }
          $i=$#coords;
          }
        }
      print "$ctg $j $score\n";
      }
    }
    @coords=();
    $ctg=$f[0];
}
      

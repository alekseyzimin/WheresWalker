#!/usr/bin/env perl
#this code computes intervals of relative homozygosity
my @diffs=();
my $ma_window=100000;
my $step_size=10000;
my $step_index=int($ma_window/ $step_size);
my $ctg="";
my $threshold=2;
my $min_window=4000000;
print "Contig Window_start Window_end\n";
while($line=<STDIN>){
chomp($line);
my @f=split(/\s+/,$line);
if(not($f[0] eq $ctg)){
  if($#diffs>0){
    detect_windows();
  }
  @diffs=();
  @coords=();
  }
  $ctg=$f[0];
  push(@diffs,$f[3]-$f[2]);
  push(@mutfreq,$f[2]);
  push(@coords,$f[1]);
}

sub detect_windows {
  my @ma=();
  for($i=$step_index;$i<$#diffs-$step_index;$i++){
    my $dave=0;
    my $mave=0;
    for($j=$i-$step_index;$j<=$i+$step_index;$j++){
      $dave+=$diffs[$j];
      $mave+=$mutfreq[$j];
    }
    $dave=$dave/($step_index*2+1);
    $mave=$mave/($step_index*2+1);
    push(@ma,$dave/(2+$mave));
  }
  #we computed MA's, now let's detect windows
  my $window_start=-1;
  my $window_end=-1;
  for($i=$step_index;$i<$#diffs-$step_index;$i++){
    if($ma[$i-$step_index]>=$threshold && $window_start==-1){
      $window_start=$coords[$i];
      $window_end=-1;
    }
    if($ma[$i-$step_index]<0 && $window_start>-1){
      $window_end=$coords[$i];
      print "$ctg $window_start $window_end\n" if($window_end-$window_start>=$min_window);
      $window_start=-1;
    }
  }
} 




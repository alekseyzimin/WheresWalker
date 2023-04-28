#!/usr/bin/env perl
#this code computes intervals of relative homozygosity
my @diffs=();
my $ma_window=750000;
my $step_size=25000;
my $step_index=int($ma_window/ $step_size);
my $ctg="";
my $threshold=0;
my $min_window=4000000;
my @intervals=();
my @lines=();
my @ma=();
my %ma_arr=();
#$threshold=$ARGV[0] if($ARGV[0]>0);
print "Contig Window_start Window_end\n";

while($line=<STDIN>){
  chomp($line);
  push(@lines,$line);
}

my $first_pass=1;
while(scalar(@intervals)==0){
  foreach $line(@lines){
    my @f=split(/\s+/,$line);
    if(not($f[0] eq $ctg)){
      if($#diffs>0){
        if($first_pass){
          @ma=();
          compute_ma();
          $ma_arr{$ctg}=[@ma];
          for($i=$step_index;$i<$#diffs-$step_index;$i++){
            print STDERR "$ctg $coords[$i] $ma[$i-$step_index]\n";
          }
        }else{
          detect_windows(); 
        }
      }
      @diffs=();
      @coords=();
    }
    $ctg=$f[0];
    push(@diffs,$f[3]-$f[2]);
    push(@mutfreq,$f[2]);
    push(@coords,$f[1]);
  }
  #print "threshold = $threshold intervals =  $#intervals\n";
  $first_pass=0;
  $threshold-=.001;
  last if($threshold < 0.02);
}

foreach my $in(@intervals){
  print "$in\n";
}

sub compute_ma {
  for($i=$step_index;$i<$#diffs-$step_index;$i++){
    my $dave=0;
    my $mave=0;
    for($j=$i-$step_index;$j<=$i+$step_index;$j++){
      $dave+=$diffs[$j];
      $mave+=$mutfreq[$j];
    }
    $dave=$dave/($step_index*2+1);
    $mave=$mave/($step_index*2+1);
    my $ma_value=$dave/(2+$mave);
    $threshold=$ma_value if($ma_value>$threshold);
    push(@ma,$ma_value);
  }
}

sub detect_windows {
  my $window_start=-1;
  my $window_end=-1;
  my @ma_local=@{$ma_arr{$ctg}};
  my $min_value=0.01;
  for($i=$step_index;$i<$#diffs-$step_index;$i++){
    if($ma_local[$i-$step_index]>=$threshold && $window_start==-1){
      my $j=$i;
      while($ma_local[$j-$step_index]>$min_value && $j>=$step_index){
        $j--;
      }
      $window_start=$coords[$j];
      $window_end=-1;
    }
    if($ma_local[$i-$step_index]<=$min_value && $window_start>-1){
      $window_end=$coords[$i];
      $window_start=$window_start-$ma_window;
      $window_start=0 if($window_start<0);
      $window_end=$window_end+$ma_window;
      #print "candidate window $ctg $window_start $window_end\n";
      push(@intervals,"$ctg $window_start $window_end") if($window_end-$window_start>=$min_window+2*$ma_window);
      $window_start=-1;
    }
  }
  #check if window started, but did not end
  if($window_start > 0){
    $window_end=$coords[-1];
    #print "end candidate window $ctg $window_start $window_end\n";
    push(@intervals,"$ctg $window_start $window_end") if($window_end-$window_start>=$min_window+2*$ma_window);
  }
} 




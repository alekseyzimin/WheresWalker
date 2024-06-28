#!/usr/bin/env perl
#this code computes intervals of relative homozygosity
my @diffs=();
my $ma_window=750000;
my $step_size=10000;
my $step_index=int($ma_window/ $step_size);
my $ctg="";
my $threshold=0;
my $min_window=$ARGV[0];
my @intervals=();
my @lines=();
my @ma=();
my %ma_arr=();
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
      @mutfreq=();
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
  last if($threshold < 0.01);
}

my $max_interval=0;
my $max_interval_index=-1;
for(my $i=0;$i<=$#intervals_sizes;$i++){
  if($intervals_sizes[$i]>$max_interval){
    $max_interval=$intervals_sizes[$i];
    $max_interval_index=$i;
  }
  #output the biggest interval
  print "$intervals[$max_interval_index]\n";
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
    #print STDERR "DEBUG $ctg $dave $mave $ma_value\n";
    $threshold=$ma_value if($ma_value>$threshold);
    push(@ma,$ma_value);
  }
}

sub detect_windows {
  my $window_start=-1;
  my $window_end=-1;
  my @ma_local=@{$ma_arr{$ctg}};
  my $min_value=0.005;
  #print "DEBUG $ctg $threshold\n";
  for($i=$step_index;$i<$#diffs-$step_index;$i++){
    if($ma_local[$i-$step_index]>=$threshold && $window_start==-1){
      my $j=$i;
      #print "DEBUG found start $ctg $threshold $coords[$i]\n";
      while($ma_local[$j-$step_index]>$min_value && $j>=$step_index){
        $j--;
      }
      #print "DEBUG went back $ctg $threshold $coords[$j]\n";
      $window_start=$coords[$j];
      $window_end=-1;
    }
    if($ma_local[$i-$step_index]<=$min_value && $window_start>-1){
      $window_end=$coords[$i];
      $window_start=$window_start-$ma_window;
      $window_start=0 if($window_start<0);
      $window_end=$window_end+$ma_window;
      #print "candidate window $ctg $window_start $window_end $min_window $ma_window\n";
      push(@intervals,"$ctg $window_start $window_end") if($window_end-$window_start>=$min_window+2*$ma_window);
      push(@intervals_sizes,$window_end-$window_start) if($window_end-$window_start>=$min_window+2*$ma_window);
      $window_start=-1;
    }
  }
  #check if window started, but did not end
  if($window_start > 0){
    $window_end=$coords[-1];
    #print "end candidate window $ctg $window_start $window_end\n";
    push(@intervals,"$ctg $window_start $window_end") if($window_end-$window_start>=$min_window+2*$ma_window);
    push(@intervals_sizes,$window_end-$window_start) if($window_end-$window_start>=$min_window+2*$ma_window);
  }
} 




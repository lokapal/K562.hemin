#!/usr/bin/perl 

my @chrs;
my @coords;

open(FILE,$ARGV[0]);
my $i=0;
while(<FILE>) {
      my @words=split('\t');
      $chrs[$i]=$words[2];
      $coords[$i]=$words[3];
      $i++;
}
close(FILE);

$start=0;

open(FILE,$ARGV[1]);
$i=0;
while(<FILE>) {
      chomp;
      chomp;
      my @words=split('\t',$_);

checkagain:   if(!(($words[0] eq $chrs[$i]) && ($coords[$i]>=$words[1]) && ($coords[$i]<=$words[2])))  {
          $i++; 
          goto checkagain;                                                                             }
      my $counts=0;
nextif:  if(($words[0] eq $chrs[$i]) && ($coords[$i]>=$words[1]) && ($coords[$i]<=$words[2]))  {
            $counts++; 
            $i++;
            goto nextif;                                                                       }
      print $_,"\t",$counts,"\n";
}
close(FILE);

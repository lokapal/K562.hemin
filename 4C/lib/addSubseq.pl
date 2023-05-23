#!/usr/bin/perl 

use Bio::DB::Fasta;

  my $db      = Bio::DB::Fasta->new($ARGV[0]);
  my @words =();
  open(FILE,$ARGV[1]);
  while(<FILE>){
        chop($_);
        @words=split(/\t/,$_);
        my $seq     = $db->seq($words[0],$words[1] => $words[2]);
        print $_,"\t",$seq,"\n";
  }
  close(FILE);

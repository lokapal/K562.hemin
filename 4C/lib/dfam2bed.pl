#!/usr/bin/perl

use strict;
use PerlIO::gzip;

#columns from DFAM that should be converted to BED
my $chrF=1;
my $begF=10;
my $endF=11;
my $strF=9;

   my $ARGC=scalar @ARGV;

   if ($ARGC==0||$ARGC>1) { die "Usage: dfam2bed dfambase.hits.gz\n"; }

    my $infile=$ARGV[0];
    open (INP, "<:gzip(autopop)", $infile) || die "Can't open input DFAM file \"$infile\": $!"; 

    my $basename=substr($infile,0,rindex($infile,"."));
    my $outfile=$basename.".bed";
    open (OUTP,">$outfile") || die "Can't create \"$outfile\": $!";

    $chrF--;
    $begF--;
    $endF--;
    $strF--;

    while(<INP>) {
       chomp;
       chomp;
       if ($_=~ m/^(\s+)?#/) { next; }  # commented strings eliminated
       if (length($_)<5) { next; }
       my @arr=split('\s+');

       if ($arr[$begF] > $arr[$endF]) { my $tmp = $arr[$begF];
                                        $arr[$begF] = $arr[$endF];
                                        $arr[$endF] = $tmp; }
       if ($arr[$chrF] =~ m/_|random/) { next; }
#       $arr[$chrF] =~s/^chr//;
       print OUTP "$arr[$chrF]\t$arr[$begF]\t$arr[$endF]\t$arr[$strF]\n";
                 }
    close (INP);
    close (OUTP);

    my $cmd="sort -k1,1V -k2,2n $outfile -o $outfile";
   `$cmd`;
   my $outfinal = "merged.bed";
   $cmd="bedtools merge -d 1 -i $outfile > $outfinal";
   `$cmd`;
   $cmd="mv $outfinal $outfile";
   `$cmd`;

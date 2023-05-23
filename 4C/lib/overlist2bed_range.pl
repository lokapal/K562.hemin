#!/usr/bin/perl
#name of file for output
    use strict;
    use POSIX;
    use List::Util qw(min max);

    my $ARGC=scalar @ARGV;

    if ($ARGC==0||$ARGC>1) { die "Usage: compareFT.pl infile.txt \n"; }

    my $infile=$ARGV[0];

    my $basename=substr($infile,0,rindex($infile,"."));

    my $outfile=$basename . ".bed";

    open (EP,$infile) || die "can't open \"$infile\" for reading: $!";
    open (OUTP,">$outfile") || die "Can't create \"$outfile\": $!";

    while (<EP>) {
         chomp;
         chomp;
         if (length($_)<5) { next; }
#chr1	21702	21755	1	chr1	21702	21755	3
         my @arr=split('\t',$_);
         my (@starts, @ends);
         my $chrom=$arr[0];
         push (@starts, $arr[1], $arr[5]);
         push (@ends, $arr[2], $arr[6]);
         my $start=min(@starts)-1;
         my $end=max(@ends)+1;
         print OUTP "$chrom\t$start\t$end\n";
         undef @starts;
         undef @ends;
                 }
    close (EP);
    close (OUTP);

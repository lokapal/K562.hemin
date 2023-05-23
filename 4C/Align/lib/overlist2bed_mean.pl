#!/usr/bin/perl
    use strict;
    use POSIX;
    use List::Util qw(min max);

    my $ARGC=scalar @ARGV;

    if ($ARGC==0||$ARGC>1) { die "Usage: overlist2bed.pl infile.txt \n"; }

    my $infile=$ARGV[0];

    my $basename=substr($infile,0,rindex($infile,"."));

    my $outfile=$basename . ".bed";
    my $outfile2=$basename . ".bedGraph";

    open (EP,$infile) || die "can't open \"$infile\" for reading: $!";
    open (OUTP,">$outfile") || die "Can't create \"$outfile\": $!";
    open (OUT2,">$outfile2") || die "Can't create \"$outfile2\": $!";

    while (<EP>) {
         chomp;
         chomp;
         if (length($_)<5) { next; }
         my @arr=split('\t',$_);
#chr1	579389	579413	5	chr1	579389	579413	3
         my (@starts, @ends);
         my $chrom=$arr[0];
         push (@starts, $arr[1], $arr[5]);
         push (@ends, $arr[2], $arr[6]);
         my $start=min(@starts);
         my $end=max(@ends);
         my $eff=($arr[3]+$arr[7])/2.0; # mean value
         print OUTP "$chrom\t$start\t$end\t$eff\n";
         print OUT2 "$arr[0]\t$arr[1]\t$arr[2]\t$arr[3]\n$arr[4]\t$arr[5]\t$arr[6]\t$arr[7]\n";
         undef @starts;
         undef @ends;
                 }
    close (EP);
    close (OUTP);
    close (OUT2);

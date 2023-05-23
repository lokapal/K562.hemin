#!/usr/bin/perl
    use POSIX;

#Column 5 = coverage
#Column 6 = reads
$COLUMN=6;

    $ARGC=scalar @ARGV;

    if ($ARGC<1) { print "Usage: table2bedg tablefile.txt [threshold]\n";
                   print "This program will convert to bedfile tablefile from alignments. Input file format:\n";
                   print "chrname<TAB>begin<TAB>end<TAB>length<TAB>coverage<TAB>reads<TAB>sequence\n";
                   exit (-1);
                  }

    $infile=$ARGV[0];

    open (INP,"$infile") || die "Can't read \"$infile\": $!";
    $basename=substr($infile,0,rindex($infile,"."));
    $outfile=$basename.".bedGraph";

    open (OUTP,">$outfile") || die "Can't create \"$outfile\": $!";

    my $threshold=0;
    if ( ($ARGC>=2) && ($ARGV[1]>0) ) { $threshold=$ARGV[1]; }
    $COLUMN--;

    if ( ($ARGC>=3) && ($ARGV[2]>4) ) { $COLUMN=$ARGV[2]-1;  }

    if ( $ARGC>=2 ) { 
        my $effcol=$COLUMN+1;
        print "Lower threshold for column $effcol set to $threshold\n"; 
                    }

    while (<INP>) {
       chomp;
       chomp;
       my @arr=split('\t',$_);
       if ($arr[$COLUMN]<$threshold) { next; }
       print OUTP "$arr[0]\t$arr[1]\t$arr[2]\t$arr[$COLUMN]\n";
                  }
   close (INP);
   close (OUTP);
   my $cmd="bedSort $outfile $outfile";
   `$cmd`;

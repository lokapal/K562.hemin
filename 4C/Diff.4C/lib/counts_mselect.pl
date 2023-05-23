#!/usr/bin/perl
use strict;

my $infile="counts.txt";
my $outfile="counts.selected";
my $geneanno="/usr/local/genomes/hg38.gtf";
my $THRESHOLD=100;

    open ANNO, $geneanno or die $!;

    my %gene2id;
    my %id2gene;

    while (<ANNO>) {
       chomp;
       chomp;
       if ($_=~ m/^(\s+)?#/) { next; }  # commented strings eliminated
       my @arr=split('\t');
       my $ptag=$arr[2];
       my $chrname=$arr[0];
       if ($chrname =~ m/GL|MT/) { next; }
       if ($ptag ne "gene") { next; }
       my @features=split(';',$arr[8]);
       my ($geneid, $genename);
       foreach my $feature (@features) {
            my @pair=split(' ',$feature);
            my $fname=$pair[0];
            my $fcont=$pair[1];
            $fcont =~ s/\"//g;
            if ($fname eq "gene_id")     { $geneid=$fcont; }
            if ($fname eq "gene_name")   { $genename=$fcont; }
                                       }
       $gene2id{$genename}=$geneid;
       $id2gene{$geneid}=$genename;
               }
    close (ANNO);


    my @files = <*.list>;
    my %ids_scored;
    foreach my $infile (@files)  { 
        open (INP,"$infile") || die "Can't read \"$infile\": $!";
        while (<INP>) {
            chomp;
            chomp;
            if ($_=~ m/^(\s+)?#/) { next; }  # commented strings eliminated
#mir-10404-1,4787171.5
            my @arr=split(',');
            my $genereads=$arr[1];
            my $geneid=$gene2id{$arr[0]};
            if (!defined $ids_scored{$geneid}) { $ids_scored{$geneid}=$arr[1]; }
                else { if ($ids_scored{$geneid} < $arr[1]) { $ids_scored{$geneid} = $arr[1]; } }
                      }
        close (INP);
                                }

    open (INP, $infile) || die "Can't open input file \"$infile\": $!"; 
    open (OUTP, ">$outfile") || die "Cannot write to file $outfile!\n";
    while (<INP>) {
       chomp;
       chomp;
       if ($_=~ m/^(\s+)?#/) { print OUTP "$_\n"; next; }  # commented strings eliminated
       my @arr=split('\t');
       my $geneid=$arr[0];
       if ($geneid eq "Geneid") { print OUTP "$_\n"; next; }
       if (defined $ids_scored{$geneid}) {
            if ($ids_scored{$geneid} >= $THRESHOLD) { 
                    print OUTP "$_\n"; 
                                                    }
                                          }
                   }
     close (INP);
     close (OUTP);

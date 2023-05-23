#!/usr/bin/perl
# this utility ignores chrM and other incomplete chromosomes with "-" or "_" in name

    use strict;
    use Fcntl qw(SEEK_SET SEEK_CUR SEEK_END);
    use File::Copy;
    use Math::Round qw(:all);

    my %exclusions = (
#         "RNA5-8SN1"  => 1, "ENSG00000278189" => 1,
#         "RNA5-8SN2"  => 1, "ENSG00000278233" => 1,
#         "RNA5-8SN3"  => 1, "ENSG00000275215" => 1,
#         "RNA5-8SP2"  => 1, "ENSG00000200434" => 1,
#         "RNA5-8SP6"  => 1, "ENSG00000251705" => 1,
#         "FP671120.1" => 1, "ENSG00000277437" => 1,
#         "FP671120.2" => 1, "ENSG00000278996" => 1,
#         "FP671120.4" => 1, "ENSG00000280800" => 1,
#         "FP236383.1" => 1, "ENSG00000280441" => 1,
#         "FP236383.2" => 1, "ENSG00000280614" => 1,
#         "FP236383.3" => 1, "ENSG00000281181" => 1,
#         "RF00002"    => 1, "ENSG00000283274" => 1,
#         "AC010970.1" => 1, "ENSG00000225840" => 1,
#         "FP671120.5" => 1, "ENSG00000281383" => 1
                     );

# Column with expression in TPM from Rsem made file
    my $COLUMN=2;
# Columnt with gene ID in the table file
    my $GENE_NAME=9;
    my $GENE_COL=10;

    my $ARGC=scalar @ARGV;

    if ($ARGC!=2) { die "Usage: tbl_addexpr tablefile.txt ExpressionFile_in_Rsem_format\n"; }

    my $infile=$ARGV[1];
    my $tablefile=$ARGV[0];

    open (INP,"$infile") || die "Can't read \"$infile\": $!";

    $COLUMN--;
    $GENE_NAME--;
    $GENE_COL--;

# First fill hashes with TSS
    my %exprhash;

    $_=<INP>;
    while (<INP>) {
       chomp;
       chomp;
       if (length($_) < 5) { next; }
       my @arr=split('\t');

       my $gene=$arr[0];
       my $expr=$arr[$COLUMN];

       $exprhash{$gene}=$expr;
                  }

    close (INP);

    my $outfile=$tablefile . ".expr";
    my $unifile=$tablefile . ".unique";
    my $allfile=$tablefile . ".all";
    open (OUTP,">$outfile") || die "Can't create \"$outfile\": $!";
    open (OUTU,">$unifile") || die "Can't create \"$unifile\": $!";
    open (OUTA,">$allfile") || die "Can't create \"$allfile\": $!";
    open (INP,"$tablefile") || die "Can't read \"$tablefile\": $!";
    my %present;
    while (<INP>) {
       chomp;
       chomp;
       if (length($_) < 5) { next; }
       my @arr=split('\t');

       my $id=$arr[$GENE_COL];
       my $name=$arr[$GENE_NAME];
       if (! defined $exprhash{$id}) { die "Gene $id is missed in the expression table?!\n"; }
       if (! defined $present{$id})  { $present{$id}=1; 
                                       if (! defined $exclusions{$id}) { print OUTU "$id\t$name\t$exprhash{$id}\n"; }
                                     }
       if (! defined $exclusions{$id}) { print OUTA "$id\t$name\t$exprhash{$id}\n"; }
       print OUTP "$_\t$exprhash{$id}\n";
                  }
    close (INP);
    close (OUTP);
    close (OUTU);
    close (OUTA);

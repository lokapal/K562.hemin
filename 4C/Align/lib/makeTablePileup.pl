#!/usr/bin/perl


open(FILE, $ARGV[0]);

$line = <FILE>;
$prevChr = "none";
$length=0;
$prevPos = -100;
while(defined($line))
{
        @words=split(/\t/,$line);

        if(($words[0] ne $prevChr) || ($words[1] ne $prevPos+1))
        {
                if($prevChr ne "none")
                {
                        $coverage = $sumCoverage/($length);
			#if($coverage > 1)
			{ print "$chr\t$start\t$end\t$length\t$coverage\n";}
                }
                        $chr = $words[0];
                        $start = $words[1];
                        $end = $words[1];
                        $sumCoverage = $words[3];
                        $length = 1;

        }
        if(($words[0] eq  $prevChr) && ($words[1] eq $prevPos+1))
        {
	  ++$length;
          $end = $words[1];
          $sumCoverage = $sumCoverage+$words[3];
        }

        $prevChr = $words[0];
        $prevPos = $words[1];
	$line = <FILE>;
}
close(FILE);


$coverage = $sumCoverage/($length);
#if($coverage > 1)
{ print "$chr\t$start\t$end\t$length\t$coverage\n";}

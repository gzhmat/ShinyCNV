#!/usr/bin/perl

my $title = join( "\t",
	( 'Name', 'Chromosome', 'Position', 'LRR', 'BAF' ) );

my $inputFile = $ARGV[0];

die("
Sex chr should be chrX and chrY
Input file name is needed!

Input file is the Final report from GenomeStudio, as follow:

[Header]
GSGT Version	1.9.4
Processing Date	3/27/2012 3:36 PM
Content		Human660W-Quad_v1_A.bpm
Num SNPs	657366
Total SNPs	657366
Num Samples	28
Total Samples	28
File 	1 of 28
[Data]
Sample ID	SNP Name	Chr	Position	Log R Ratio	B Allele Freq
CHA1	200003	9	139026180	0.0863	0.0136
CHA1	200006	9	139046223	0.1167	0.9915
CHA1	200047	2	219793146	0.2250	0.5685
....
..

") unless ($inputFile);
my $outputFile = '';
print "Processing $inputFile!\n";
my %hs = ();
open( IN, "$inputFile" );
while (<IN>) {
	chomp;
	s/\r//g;
	my ( $sample, $snp, $chr, $poz, $lrr, $baf ) = split(/\t/);
	if ($baf=~/\d+/ && $chr ne 0 && $chr =~ /\d|[XY]$/i ) {
		$chr = $chr eq "XY" ? "X" : $chr;
		$hs{$chr}{$poz} = [ $snp, $lrr, $baf ];
		$outputFile = $sample . '.txt' unless ($outputFile);
	}
}
close IN;

print "Hash is over!\n";

open( OUT, ">$outputFile" );
print OUT $title, "\n";

for my $chr ( sort chrSort keys %hs ) {
	my @pozs = keys %{$hs{$chr}};
	@pozs = sort { $a <=> $b } @pozs;
	for my $poz (@pozs) {
		my ( $snp, $lrr, $baf ) = @{ $hs{$chr}{$poz} };
		print OUT join( "\t", ( $snp, $chr, $poz, $lrr, $baf ) ), "\n";
	}
}
close OUT;

sub chrSort {
	my ($ai,$bi)=($a,$b);
	$ai = $ai eq 'X' ? 23 : $ai eq 'Y' ? 24 : $ai;
	$bi = $bi eq 'X' ? 23 : $bi eq 'Y' ? 24 : $bi;
	$ai <=> $bi;
}

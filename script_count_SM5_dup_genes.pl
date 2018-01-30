#!/usr/bin/perl

use strict;

my $leftDupStart = 6012459;
my $leftDupEnd = 6916809;
my $rightDupStart = 21972072;
my $rightDupEnd = 22689962;
my %genes;

open INF,"dmel-all-r6.17.gtf";
while (<INF>) {
	my(@F) = split /\t/, $_;
	next unless $F[0] eq "2R";
	my($geneID) = $F[8] =~ /gene_id \"(\w+)\"\;/;
	my($geneSymbol) = $F[8] =~ /gene_symbol \"(\S+)\"\;/;

	my $inInterval = 0;

	foreach my $id ($F[3]..$F[4]) {
		$inInterval = 1 if $id >= $leftDupStart && $id <= $leftDupEnd;
		$inInterval = 1 if $id >= $rightDupStart && $id <= $rightDupEnd;
	}
	next if $inInterval == 0;

	if ($genes{$geneID}) {
		my($currLow,$currHigh) = $genes{$geneID} =~ /chr2R\,(\d+)\,(\d+)\,/;

		$currLow = $F[3] if $F[3] < $currLow;
		$currHigh = $F[4] if $F[4] > $currHigh;

		$genes{$geneID} = "chr2R,$currLow,$currHigh,$geneID,$geneSymbol";
	} else {
		$genes{$geneID} = "chr2R,$F[3],$F[4],$geneID,$geneSymbol";
	}
	
}
close INF;

my($count,$sno,$trna,$cr);
foreach (keys %genes) {
	print "$genes{$_}\n";
	++$count;
	$sno++ if $genes{$_} =~ /sno/;
	$cr++ if $genes{$_} =~ /CR\d+/;
	$trna++ if $genes{$_} =~ /trna/i;
}
my $totalGenes = $count - $sno - $cr - $trna;

print "Genes: $totalGenes\n";
print "Sno: $sno\n";
print "trna: $trna\n";
print "CR: $cr\n";

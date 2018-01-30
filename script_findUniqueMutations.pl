#!/usr/bin/perl

use strict;

my %genePosition;
open INF,"dmel-all-r6.17.gtf";
while (<INF>) {
        my(@F) = split /\t/, $_;
        next unless $F[0] eq "2R" || $F[0] eq "2L";
        my($geneSymbol) = $F[8] =~ /gene_symbol \"(\S+)\"\;/;
	$F[0] = "chr$F[0]";

        foreach my $id ($F[3]..$F[4]) {
		$genePosition{$F[0]}{$id} = $geneSymbol;
        }
}
close INF;

my(%singletons,%sGenes);
foreach my $foo ("SM5","SM6a","CyO") {
	open INF,"$foo.singletons.out";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		my(@D) = split /\//, $F[4];
		my $line = $D[6];
		$singletons{$F[0]}{$F[1]} = $line;
	}
}

my %mut;
open INF,"singletons.ann.vcf";
while (<INF>) {
	chomp($_);
	my(@F) = split /\t/, $_;
	next unless $singletons{$F[0]}{$F[1]};
	my($ann) = $F[7] =~ /ANN\=(.+)/;
	foreach my $foo (split /\,/, $ann) {
		next if $foo =~ /(upstream_gene_variant|downstream_gene_variant|intragenic_variant|intergenic_region|5_prime_UTR_variant|3_prime_UTR_variant|synonymous_variant|non_coding_transcript_exon_variant|5_prime_UTR_premature_start_codon_gain_variant|splice_region_variant&intron_variant|splice_region_variant&non_coding_transcript_exon_variant|stop_lost)/;
		next if $foo =~ /intron_variant/ && $foo !~ /splice/;
		my(@bar) = split /\|/, $foo;
                my $gene = $genePosition{$F[0]}{$F[1]};
		my $stock = $singletons{$F[0]}{$F[1]};
                #print "$singletons{$F[0]}{$F[1]}\t$F[0]\t$F[1]\t$bar[1]\t$gene\t$F[3] -> $F[4]\t$bar[10]\n";
                $mut{$F[0]}{$F[1]}{$stock} = "$F[0]\t$F[1]\t$bar[1]\t$gene\t$F[3] -> $F[4]\t$bar[10]";
		#$sGenes{$bar[3]} = "$F[0]|$F[1]";
		#print "$bar[3]\t$F[0]\t$F[1]\t$foo\n";
	}
}
close INF;

foreach my $chr (keys %mut) {
        foreach my $id (keys %{$mut{$chr}}) {
                foreach my $stock (sort keys %{$mut{$chr}{$id}}) {
                        print "$stock\t$mut{$chr}{$id}{$stock}\n";
                }
        }
}

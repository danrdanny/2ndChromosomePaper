#!/usr/bin/perl

use strict;

my $output;
foreach my $stock ("SM5","SM6a","CyO","shared","singletons") {
	my(%singletons,%sGenes);
	if ($stock =~ /singletons/) {
		foreach my $foo ("SM5","SM6a","CyO") {
			open INF,"$foo.singletons.out";
			while (<INF>) {
				chomp($_);
				my(@F) = split /\t/, $_;
				my(@D) = split /\//, $F[4];
				my $line = $D[6];
				$singletons{$F[0]}{$F[1]} = $line;
			}
			close INF;
		}

		open INF,"singletons.ann.vcf";
		while (<INF>) {
			chomp($_);
			my(@F) = split /\t/, $_;
			next unless $singletons{$F[0]}{$F[1]};
			my($ann) = $F[7] =~ /ANN\=(.+)/;
			foreach my $foo (split /\,/, $ann) {
				next if $foo =~ /(upstream_gene_variant|downstream_gene_variant|intragenic_variant|intergenic_region|5_prime_UTR_variant|3_prime_UTR_variant|synonymous_variant)/;
				next if $foo =~ /intron_variant/ && $foo !~ /splice/;
				my(@bar) = split /\|/, $foo;
				$sGenes{$bar[3]} = "$F[0]|$F[1]";
				#print "$bar[3]\t$F[0]\t$F[1]\t$foo\n";
			}
		}
		close INF;
	}

	my(%genes);
	my $file = "shared.$stock.genes.txt";
	   $file = "shared.genes.txt" if $stock =~ /shared/i;
	   $file = "singletons.genes.txt" if $stock =~ /singletons/i;
	open INF,"$file";
	while (<INF>) {
		my(@F) = split /\t/, $_;
		my $gene = $F[0];

		$genes{"prematureStart"}{$gene} = 1 if $F[9] > 0; # premature start
		$genes{"missense"}{$gene} = 1 if $F[13] > 0; #missense;
		$genes{"spliceAcceptor"}{$gene} = 1 if $F[16] > 0; #splice acceptor
		$genes{"stop"}{$gene} = 1 if $F[17] > 0; #stop
	}
	close INF;

	foreach my $class ("missense","spliceAcceptor","stop") {
		my($count,$geneList);
		foreach my $gene (keys %{$genes{$class}}) {
			$geneList .= "$gene ";
			++$count;

			$output .= "$stock\t$class\t$gene\t";

			if ($stock =~ /singletons/) {
				my($chr,$id) = $sGenes{$gene} =~ /(\w+)\|(\d+)/;
				$output .= "$singletons{$chr}{$id}\n"; 
			} else {
				$output .= "\n";
			}
		}
		print "$stock\t$class\t$count\n";
	}
}

print "$output\n";

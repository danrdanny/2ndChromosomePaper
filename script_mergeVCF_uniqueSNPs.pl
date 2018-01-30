#!/usr/bin/perl

use strict;
use Getopt::Std;

my $chrL = "chr2L";
my $chrR = "chr2R";

# Get merged data
my(%mergedData,%allData,%singletons,%allDataCount,$singletonSNPvcfData);
foreach my $stock ('CyO','SM5','SM6a') {
	my $filename = "$stock.merged.vcf.gz";
	print "[".localtime()."] opening $filename\n";
       	open(INF,sprintf("zcat < %s |", $filename)) or die "Can't open $filename: $!";
       	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		next unless $F[0] eq $chrL || $F[0] eq $chrR;
		next if $_ =~ /INDEL/;
		next unless $_ =~ /0\/1/;
		next unless $F[5] > 220;

		$mergedData{$stock}{$F[0]}{$F[1]} = $_;
		$allDataCount{$F[0]}{$F[1]}++;
		$allData{$F[0]}{$F[1]} = $_;
	}
	close INF;

	open INF,"$stock.singletons.out" or die "can't open $stock singletons: $!";
	while (<INF>) {
		chomp($_);
		my(@F) = split /\t/, $_;
		$singletons{$stock}{$F[0]}{$F[1]} = $_;
		
		$singletonSNPvcfData .= "$allData{$F[0]}{$F[1]}\n";
	}
	close INF;
}

my($totalSNPcount,$sharedSNPvcfData);
foreach my $chr (keys %allData) {
	foreach my $id (keys %{$allData{$chr}}) {
		next if $allDataCount{$chr}{$id} == 1;
		$totalSNPcount++; 

		$sharedSNPvcfData .= "$allData{$chr}{$id}\n";
	}
}
print "Total SNPs: $totalSNPcount\n";
`zcat < SM6a.merged.vcf.gz | grep "#" | grep -v "tools" > shared.vcf`;
open OUTF,">>shared.vcf";
print OUTF $sharedSNPvcfData;
close OUTF;

`zcat < SM6a.merged.vcf.gz | grep "#" | grep -v "tools" > singletons.vcf`;
open OUTF,">>singletons.vcf";
print OUTF $singletonSNPvcfData;
close OUTF;


my @SM5stocks = qw/ SM5-1143 SM5-223 SM5-240 SM5-400 SM5-405 /;
my @CyOstocks = qw/ CyO-1602 CyO-2 CyO-3076 CyO-31 CyO-471 CyO-533 CyO-TM3-504 CyO-TM3-mp-22239 SM6a-8785 /;
my @SM6astocks = qw/ SM6a-1465 SM6a-23663 SM6a-24380 SM6a-6853 SM6a-9162 /;
my @ALLstocks = qw/ SM5-1143 SM5-223 SM5-240 SM5-325 SM5-400 SM5-405 CyO-1602 CyO-2 CyO-3076 CyO-31 CyO-471 CyO-533 CyO-TM3-24759 CyO-TM3-504 CyO-TM3-mp-22239 SM6a-1465 SM6a-23663 SM6a-24380 SM6a-6853 SM6a-8785 SM6a-9162 /;

my %stockOutput;
foreach my $stock ('CyO','SM5','SM6a') {
	my($count,$sharedSNPvcfData);
	foreach my $chr (keys %{$mergedData{$stock}}) {
		foreach my $id (keys %{$mergedData{$stock}{$chr}}) {
			next unless $allDataCount{$chr}{$id} == 1;
			next if $singletons{$stock}{$chr}{$id};

			$sharedSNPvcfData .= "$allData{$chr}{$id}\n";
			#print "$stock Unique: $chr $id\n";
			++$count;
		}
	}
	print "SNPs unique to $stock: $count\n";
	#$totalSNPcount -= $count;
	`zcat < SM6a.merged.vcf.gz | grep "#" | grep -v "tools" > shared.$stock.vcf`;
	open OUTF,">>shared.$stock.vcf";
	print OUTF $sharedSNPvcfData;
	close OUTF;
}

#print "Total shared SNPs: $totalSNPcount\n";

exit 0;

my(%vcfData,%snpCount,$stockCount,@stocks,$pwd,$vcfFileCount);
foreach my $stock (@stocks) {
	if (-e "$pwd/$stock/$stock.vcf.gz") {
		++$stockCount;
		my $filename = "$stock/$stock.vcf.gz";
		print "[".localtime()."] opening $filename\n";
        	open(INF,sprintf("zcat %s |", $filename)) or die "Can't open $filename: $!";
        	while (<INF>) {
			my(@F) = split /\t/, $_;
			next unless $F[0] eq $chrL || $F[0] eq $chrR;
			next if $_ =~ /INDEL/;
			next unless $_ =~ /0\/1/;
			#next unless $F[4] =~ /^(A|G|C|T)$/;
			next unless $F[5] > 220;

			$snpCount{$F[0]}{$F[1]}++;
			$vcfData{$stock}{$F[0]}{$F[1]}++; # = $F[4];
		}
		close INF;

		++$vcfFileCount;
	} else {
		print "[".localtime(time)."] WARN: Atttempting multi-sample vcf calling and $stock.vcf is missing.\n";
	}

}

#open OUTF,">$pwd/out_mergedVCFs/$opts{'c'}.singletons.tsv";
#print OUTF $singletonOutput;
#close OUTF;


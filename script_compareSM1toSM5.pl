#!/usr/bin/perl

use strict;
use Getopt::Std;

## User-defined variables
my $minVCFScore = 200;

my $numberToReport = 1; # marks unique sites less than or equal to this number

my $chrL = "chr2L";
my $chrR = "chr2R";

#@stocks = qw/ SM6a-1465 SM6a-23663 SM6a-24380 SM6a-6853 SM6a-9162 /;
my $sm6Stock = "SM5-223";
my @stocks = qw/ SM5-325 CyO-TM3-24759 SM6a-TM3-mp-lab /;

my %chrSizes;
   $chrSizes{'chrX'} = 23542271;
   $chrSizes{'chr2L'} = 23513712;
   $chrSizes{'chr2R'} = 25286936;
   $chrSizes{'chr3L'} = 28110227;
   $chrSizes{'chr3R'} = 32079331;
   $chrSizes{'chr4'} = 1348131;

## Get our current working directory
my $pwd = `pwd`; chomp($pwd);


my(%refVCFData); 
my $filename = "$sm6Stock/$sm6Stock.vcf.gz";
print "[".localtime()."] opening $filename\n";
open(INF,sprintf("gunzip -c $filename |")) or die "Can't open $filename: $!";
while (<INF>) {
	my(@F) = split /\t/, $_;
	next unless $F[0] eq $chrL || $F[0] eq $chrR;
	next if $_ =~ /INDEL/;
	next unless $_ =~ /0\/1/;
	next unless $F[4] =~ /^(A|G|C|T)$/;

	#next unless $F[5] > 120;
	$refVCFData{$F[0]}{$F[1]} = $F[4];
}
close INF;

my(%vcfData);
foreach my $stock (@stocks) {
	my $filename = "$stock/$stock.vcf.gz";
	print "[".localtime()."] opening $filename\n";
       	open(INF,sprintf("gunzip -c $filename |")) or die "Can't open $filename: $!";
       	while (<INF>) {
		my(@F) = split /\t/, $_;
		next unless $F[0] eq $chrL || $F[0] eq $chrR;
		next if $_ =~ /INDEL/;
		next unless $_ =~ /0\/1/;
		next unless $F[4] =~ /^(A|G|C|T)$/;

		next unless $F[5] > 220;
		$vcfData{$stock}{$F[0]}{$F[1]} = $F[4];
	}
	close INF;
}

my $output = "stock\tchr\tpos\tpercent\n";
foreach my $stock (@stocks) {
	print "$stock\t";
	foreach my $chr ($chrL,$chrR) {
		my $maxPos = $chrSizes{$chr};
		my $currMin 	= 0;
		my $step 	= 10000; # in kb
		my $currMax 	= 10000;
		while ($currMax < $maxPos) {
			my $count = 0;
			foreach my $id ($currMin..$currMax) { 
				next if $refVCFData{$chr}{$id} eq $vcfData{$stock}{$chr}{$id};

				if ($refVCFData{$chr}{$id} && $vcfData{$stock}{$chr}{$id} && ($refVCFData{$chr}{$id} ne $vcfData{$stock}{$chr}{$id})) {
					++$count;
				} elsif ($refVCFData{$chr}{$id}  && !$vcfData{$stock}{$chr}{$id}) {
					#++$count;
				} elsif (!$refVCFData{$chr}{$id} && $vcfData{$stock}{$chr}{$id}) {
					++$count;

					print "$stock\t$chr\t$id\t$refVCFData{$chr}{$id}\t$vcfData{$stock}{$chr}{$id}\n";
				}
			}

			$count = 0 if $count <= 1; # 1 SNP/interval gets a little noisy
			$count = 1 if $count > 1;
			#my $frac = sprintf("%0.0f",(($count / $step) * 10000));
			#$count = 0 if $frac < 0.5;
			#$count = 1 if $frac >= 0.5;
			#$count = $frac;
			#print "$stock\t$chr\t$currMax\t$count\n";

			$output .= "$stock\t$chr\t$currMax\t$count\n";

			$currMin = $currMax + 1;
			$currMax += $step;
		}
	}
}

open OUTF,">./out_SM1_vs_SM5.diff.tsv";
print OUTF $output;
close OUTF;

#!/usr/bin/perl

use strict;
use Getopt::Std;

## User-defined variables
my $minVCFScore = 200;

## Command-line options
my %opts;
getopts('c:h', \%opts); # Values in %opts

my $numberToReport = 1; # marks unique sites less than or equal to this number

my(@stocks,$chrL,$chrR);

$chrL = "chr2L";
$chrR = "chr2R";
# 22 total stocks
if ($opts{'c'} =~ /sm5/i) {

	@stocks = qw/ SM5-1143 SM5-223 SM5-240 SM5-400 SM5-405 /;

} elsif ($opts{'c'} =~ /cyo/i) {

	@stocks = qw/ CyO-1602 CyO-2 CyO-3076 CyO-31 CyO-471 CyO-533 CyO-TM3-504 CyO-TM3-mp-22239 SM6a-8785 /;

} elsif ($opts{'c'} =~ /sm6a/i) {

	@stocks = qw/ SM6a-1465 SM6a-23663 SM6a-24380 SM6a-6853 SM6a-9162 /;

} elsif ($opts{'c'} =~ /all/i) {

	@stocks = qw/ SM5-1143 SM5-223 SM5-240 SM5-325 SM5-400 SM5-405 CyO-1602 CyO-2 CyO-3076 CyO-31 CyO-471 CyO-533 CyO-TM3-24759 CyO-TM3-504 CyO-TM3-mp-22239 SM6a-1465 SM6a-23663 SM6a-24380 SM6a-6853 SM6a-8785 SM6a-9162 /;

} elsif ($opts{'c'} =~ /sm1/i) {

	@stocks = qw/ SM5-325 CyO-TM3-24759 SM6a-TM3-mp-lab /;

}

## Usage Statement
if ($opts{'h'}) {
	print "

	-c	stock
	-h	This helpful help.

	\n";
	exit 0;
}

## Subroutines
sub executeCommand {
        #open LOGF, ">>$pwd/out_log/$logfile";
        #print LOGF "[".localtime()."] CMD: $_[0]\n";
        #close LOGF;
        my $output = `$_[0]`;
        return($output);
}

my %chrSizes;
   $chrSizes{'chrX'} = 23542271;
   $chrSizes{'chr2L'} = 23513712;
   $chrSizes{'chr2R'} = 25286936;
   $chrSizes{'chr3L'} = 28110227;
   $chrSizes{'chr3R'} = 32079331;
   $chrSizes{'chr4'} = 1348131;

## Get our current working directory
my $pwd = `pwd`; chomp($pwd);

# Create common and unique VCF files for all chromosomes sequenced in this study
print "[".localtime(time)."] RUN: Creating common VCF files.\n";

my($VCFFiles,$vcfFileCount);

my(%vcfData,%snpCount,$stockCount);
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
			next unless $F[4] =~ /^(A|G|C|T)$/;

			$snpCount{$F[0]}{$F[1]}++;

			next unless $F[5] > 200;
			$vcfData{$stock}{$F[0]}{$F[1]} = $F[4];
		}
		close INF;

		#$VCFFiles .= "$pwd/out_mergedVCFs/$stock.vcf.gz ";
		#executeCommand("cat $pwd/$stock/$stock.vcf | grep \"^#\" > $pwd/out_mergedVCFs/$stock.vcf");
		#executeCommand("awk -v x=200 '\$6 > x' $pwd/$stock/$stock.vcf | grep \"^$chrL\" | grep \"0/1\" | grep -v INDEL >> $pwd/out_mergedVCFs/$stock.vcf");
		#executeCommand("awk -v x=200 '\$6 > x' $pwd/$stock/$stock.vcf | grep \"^$chrR\" | grep \"0/1\" | grep -v INDEL >> $pwd/out_mergedVCFs/$stock.vcf");
		#executeCommand("bgzip -f $pwd/out_mergedVCFs/$stock.vcf");
		#executeCommand("tabix -f -p vcf $pwd/out_mergedVCFs/$stock.vcf.gz");
		++$vcfFileCount;
	} else {
		print "[".localtime(time)."] WARN: Atttempting multi-sample vcf calling and $stock.vcf is missing.\n";
	}

}

my $commonSites = $stockCount - 1;

my(%singletons);
my $singletonOutput = "Stock\tChr\tPos\n";
foreach my $stock (@stocks) {
	print "$stock\t";
	foreach my $chr ($chrL,$chrR) {
		my $singletonCount = 0;
		foreach my $id (sort {$a<=>$b} keys %{$vcfData{$stock}{$chr}}) {
			next unless $snpCount{$chr}{$id} <= $numberToReport;
			#next if $snpCount{$chr}{$id} == $commonSites;

			++$singletonCount;
			$singletons{$stock}{$chr}{$id} = 1;
			#print "$stock\t$chr\t$id\n";
			$singletonOutput .= "$stock\t$chr\t$id\n";
		}
		print "$singletonCount\t";
	}
	print "\n";
}
open OUTF,">$pwd/out_mergedVCFs/$opts{'c'}.singletons.tsv";
print OUTF $singletonOutput;
close OUTF;

#print "[".localtime(time)."] RUN:  Merging $vcfFileCount VCF files.\n";
#executeCommand("vcf-merge $VCFFiles > $pwd/out_mergedVCFs/$opts{'c'}.merged.vcf");
#executeCommand("vcftools --singletons --vcf $pwd/out_mergedVCFs/$opts{'c'}.merged.vcf");
#executeCommand("mv $pwd/out.singletons $pwd/out_mergedVCFs/$opts{'c'}.singletons.out");
#my $singletons = executeCommand("cat $pwd/out_mergedVCFs/$opts{'c'}.singletons.out | wc -l");
#my $singletonCount = executeCommand("cut -f 5 $pwd/out_mergedVCFs/$opts{'c'}.singletons.out | sort | uniq -c");
#print "[".localtime(time)."] Singletons: $singletons";
#print "[".localtime(time)."] Per-stock Singletons:\n$singletonCount";

#---------------------------------------------------------#
#        Create singleton heatmap feeder file
#---------------------------------------------------------#

print "[".localtime(time)."] FIG: Iterating over singletons to make heatmap source file.\n";
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
			foreach my $id (sort {$a<=>$b} keys %{$singletons{$stock}{$chr}}) {
				# this shrinks the hash as we go through it, speeding the program up
				$singletons{$stock}{$chr}{$id} = undef if $id < $currMin;

				next unless $id >= $currMin && $id <= $currMax;
				++$count;
			}
			$count = 0 if $count <= 1; # 1 SNP/interval gets a little noisy
			#print "$stock\t$chr\t$currMax\t$count\n" if $count > 0;
			$count = 1 if $count > 1;

			$output .= "$stock\t$chr\t$currMax\t$count\n";

			$currMin = $currMax + 1;
			$currMax += $step;
		}
	}
}
open OUTF,">$pwd/out_mergedVCFs/$opts{'c'}.singletons.forheatmap.tsv";
print OUTF $output;
close OUTF;

#!/usr/bin/perl

## This script merges vcf files and counts both shared and unique SNPs for 
## all balancers described in "The molecular and genetic characterization of 
## second chromosome balancers in Drosophila melanogaster". 

use strict;
use Getopt::Std;

my %opts;
getopts('c:h', \%opts); # Values in %opts
my @stocks;

if ($opts{'h'} || !$opts{'c'}) {
	print "
	Required:
	-c	The chromosome you'd like to merge and do counts on, valid inputs are:
		 SM5 CyO SM6a or all

	Optional:
	none

	\n";

	exit 0;
}

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

print "[".localtime(time)."] Option $opts{'c'} passed\n";

my($gzVCFFiles,$vcfFileCount);
foreach my $stock (@stocks) {
	print "[".localtime(time)."] getting SNPs from $stock\n";
	if (-e "$stock/$stock.vcf.gz") {
		$gzVCFFiles .= "$stock.vcf.gz ";
		`zcat < $stock/$stock.vcf.gz | grep \"^#\" > $stock.vcf`;
		`zcat < $stock/$stock.vcf.gz | awk -v x=220 '\$6 > x' | grep \"0/1\" | grep -v INDEL >> $stock.vcf`;
		`bgzip -f $stock.vcf`;
		`tabix -f -p vcf $stock.vcf.gz`;
		++$vcfFileCount;
	} else {
		die "Attempting multi-sample vcf calling and $stock.vcf is missing.\n";
	}
}
		
print "[".localtime(time)."] Merging $vcfFileCount VCF files.\n";
`vcf-merge $gzVCFFiles > $opts{'c'}.merged.vcf`;
`vcftools --singletons --vcf $opts{'c'}.merged.vcf`;
`cat out.singletons | grep "chr2" > $opts{'c'}.singletons.out`;
`rm -f out.singletons`; 
my $singletons = `cat $opts{'c'}.singletons.out | wc -l`;
my $singletonCount = `cut -f 5 $opts{'c'}.singletons.out | sort | uniq -c`;
print "[".localtime(time)."] $opts{'c'} singletons: $singletons";
print "[".localtime(time)."] $opts{'c'} Singleton details:\n$singletonCount";


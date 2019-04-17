#!/usr/bin/env perl

# 	Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk)
# 	Version 0.1
#	Part of interproscan and blast2go pipeline

use strict;
use warnings;

my $usage = "\nUsage: perl $0 <input file>\n
Note: The file needed to be sorted by first column (sort -k1,1 <filename>, will do)

";
my $file = shift or die $usage;

my $inf;
my $true=1;
my $for_pass="";
my ($id,$goec,$desc);
my (@go,@ec,@desc_a);
open $inf, "<" . $file;

while (<$inf>) { #reading file and storing contig,go-ec numbers and annotation to array
	chomp();
	($id,$goec,$desc)= split(/\t/);
	if($true) {
		$for_pass=$id;
		$true=0;
	}
	if ($id ne $for_pass) {
		local $"="|";
		if (scalar(@ec) == 0) {
			print "$for_pass\t@go\tNULL\t@desc_a\n";
		}
		else {
			print "$for_pass\t@go\t@ec\t@desc_a\n";
		}
		$for_pass=$id;
		#empty array
		(@go,@ec,@desc_a)=();
	}
	if($id eq $for_pass) {
		if ($goec =~ /^GO:/) {
			push (@go,$goec);
		}
		if ($goec =~ /^EC:/) {
			push (@ec,$goec);
		}
		if ($desc) {
			push (@desc_a,$desc);
		}
	}
	if(eof($inf)) {
		local $"="|";
		if (scalar(@ec) == 0) {
			print "$for_pass\t@go\tNULL\t@desc_a\n";
		}
		else {
			print "$for_pass\t@go\t@ec\t@desc_a\n";
		}
	}
}
close($inf);

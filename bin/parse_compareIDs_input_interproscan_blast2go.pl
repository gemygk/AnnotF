#!/usr/bin/env perl

# 	Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk)
# 	Version 0.1
#	Part of interproscan and blast2go pipeline

use strict;
use warnings;

my $usage = "\n\tUsage: perl $0 <complete_contig_list[FILE]> <final.txt.out[FILE]>\n";

my $id_file = shift or die $usage;
my $gff_file = shift or die $usage;

open(ID, $id_file) || die "Cannot open $id_file\n";
open(GFF, $gff_file) || die "Cannot open $gff_file\n";
my %hash1=();

#first file --- the file having the genes with functional annotation
while (my $line=<ID>) {
	chomp($line);
	my ($contig,$agi) = split("\t",$line);
	$hash1{$contig} = $agi;
}
close(ID);

#second file --- genes for which FA is needed
#print "ID\tBlast2GO_GO_term\tBlast2GO_EC_Number\tBlast2GO_GO_Description\tInterproscan_IPR\tInterproscan_IPR_Description\tInterproscan_GO_term\tInterproscan_EC_Number\tInterproscan_Pathways\n";
while (<GFF>) {
	chomp();
	my ($first,$rest) = split("\t",$_,2);
	if (exists $hash1{$first}) { #comp0_c0        comp0_c0_seq1
		print "$first\t$rest\n";
		delete $hash1{$first};
	}
	
}
close (GFF);
foreach (keys %hash1) {
	print "$_\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\tNULL\n";
}
exit;

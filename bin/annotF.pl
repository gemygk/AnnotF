#!/usr/bin/env perl

#	This pipeline with execute interproscan, blastp for blast2go for input "peptides"
#	at the same time and finally combine all the results as a functional annotation file

#	This script use packages : LEAFF, INTERPRO, BLAST2GO (including blastp)

# 	Version 1.01
#	- Added XML,GFF3 interproscan output formats
#	- Changed the blast database location from this /tgac/references/databases/blast/nr to this /tgac/public/databases/blast/ncbi/nr

# 	Version 1.0
#	- First release

# AUTHOR: Gemy George Kaithakottil (Gemy.Kaithakottil@tgac.ac.uk || gemygk@gmail.com)


use strict;
use warnings;
use File::Basename;
use Parallel::ForkManager;
use threads;
#use Getopt::Long qw(:config no_ignore_case pass_through); 
use Getopt::Long; 
use Cwd;

my $help_flag;
my $input;
my $num = 1000;
my $queue = "Test128";
my $help;
my $prog = basename($0);

my $usage = <<_EOUSAGE_;

#########################################################################################################
#
#    annotF (v1.01) is a functional annotation pipeline for the annotation of query proteins
#   
#        Usage: $prog [options] --fasta <protein.fasta>
#
#	Required:
#	--fasta <string>		fasta file properly formatted 
#					(without any "dots" in sequences or "empty/blank lines")
#
#	General Options:
#	--chunk_size <int>		required chunk size (fasta per file) for faster execution of jobs
#					ideal chunk size would be 1000
#					(default: $num fasta sequences per file)
#	--queue <string>		queue to submit jobs (default: Test128)
#	--help				print this option menu and quit
#							
#    NOTE: 
#    Only fasta header should be present for fasta sequence and no pipe ("|") allowed in fasta header
#
#    Example of a peptide fasta file:
#    >A2YIW7
#    MAAEEGVVIACHNKDEFDAQMTKAKEAGKVVIIDFTASWCGPCRFIAPVFAEYAKKFPGAVFLKVDVDELKEVAEKYNVE
#    AMPTFLFIKDGAEADKVVGARKDDLQNTIVKHVGATAASASA
#    >comp21620_c0_seq1
#    MMGKGKVLVYHFLPSIFLVNITYYIYIYIYIYIFLYIHTYNYIHIYILLMASITGSAVSI
#    SSFSCSFKLNQASARVSTLNSVPFSINGKSFPSIKLRPAQRFQVSCMAAKPETVEKVCGI
#    VRKQLAIAADTEITGESKFAALGADSLDTVEIVMGLEEEFGISVEEESAQTIATVQDAAD
#    LIEKLLA
#
# Contact : Gemy George Kaithakottil (gemy.kaithakottil\@tgac.ac.uk || gemygk\@gmail.com)
#########################################################################################################
_EOUSAGE_
    ;

unless (@ARGV) {
    die "$usage\n";
}

&GetOptions ( 'fasta=s' => \$input,
              'chunk_size=i' => \$num,
              'queue=s' => \$queue,
              'h|help' => \$help,

);

#if ($help_flag || !$input || !$num) {
if ($help) {
    die $usage;
}

if (@ARGV) {
    die "Error, do not understand options: @ARGV\n";
}

unless ($input) {
    die "ERROR: No input fasta provided.\n$usage\n";
}

unless(-e $input) {
		print "ERROR: Cannot open $input\n$!\n";
		exit;
	}

# v1.01 update
# Link the protein file to the current location
my $input_basename = basename($input);
qx(ln -s $input .) unless (-e $input_basename);
$input = $input_basename;

#unless ($queue) {
#    die "$usage\n";
#}

	
#my $log;
#if (-e "$input.log") { print "Removing the log file already present $input.log\n";`rm $input.log`; }
#open $log, ">>" . "$input.log";

# Get working directory location
chomp(my $pwd=`pwd`);

#To pass number of chunks created in leaff step;
my $chunk;
my $prefix;

# CONSTANTS
# Number of jobs to be executed at a time -- default 500 jobs
my $THREADS_NUM = 1000; # http://bickson.blogspot.co.uk/2011/12/multicore-parser-part-2-parallel-perl.html

# Print version status:
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
print "%%          annotF - v1.01 (Release:05-Nov-2015)        %%\n";
print "%%   Gemy.Kaithakottil\@tgac.ac.uk || gemygk\@gmail.com   %%\n";
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

# Create directories
system("mkdir -p leaff interproscan blastpForBlast2go blast2go temp_folder");

#########################
##LEAFF##################
#########################
sub leaff() {
	#print $log "Inside leaff\n";
	my $file = $input;
	$prefix = "$input.pre";  # symbol $prefix used for blastp step
	my $count=0;
	my $lines;
	$lines = `grep -c '^>' $file`; # MAIN
	chomp($lines);
	my $total = ($lines/$num);
	my @commands=();
	$num = $num-1;
	for (my $i=0;$i<=($lines-1);$i++) {
		$count++;
		my $j=$i+$num; # like 9,19,29...
		if ($j <= $lines && ($j != $lines)) {
		#if ($j <= $lines) {
			#print $log "Starting ",($i+1),"-",($j+1)," End\n";
			push (@commands, "bsub -q $queue -K -o $pwd/leaff/out_leaff_$count -J leaff_$count \"source leaff-13_12_2012;leaff -f $file -S $i $j > $pwd/leaff/$prefix-$count.txt\"");
			$i=$j;
		}
		else {
			$j = ($lines-1);
			#print $log "Starting ",($i+1),"-",($j+1)," End**Last File**\n";
			push (@commands, "bsub -q $queue -K -o $pwd/leaff/out_leaff_$count -J leaff_$count \"source leaff-13_12_2012;leaff -f $file -S $i $j > $pwd/leaff/$prefix-$count.txt\"");
			last;
		}
	}
	$chunk = $count; # see whether its used again or else delete it
	
	#http://stackoverflow.com/questions/3991143/how-can-i-pass-two-arrays-and-a-string-to-a-perl-subroutine
	#pass array as reference
	&process_cmd(\@commands,"LEAFF");
	#print $log "leaff process completed!!\n";
	#print "Chunking fasta executed!!\n";
} #end of subroutine

#########################
###Interproscan_rc2######
#########################

sub interproscanRC2() {
	chomp(my @files=`ls $pwd/leaff/*.txt`); # read from leaff directory
	my $count = 1;
	my @commands=();
	foreach (@files) {
		my @input=split (/\//);
		# Actual way to run the script
		push(@commands,"bsub -q $queue -K -J iprscan5_$count -o $pwd/interproscan/out_interproscan_rc2_$input[$#input]-$count -n 8 -R \"span[ptile=8] rusage[mem=4096]\" \"source interproscan-5;interproscan.sh -i $_ -b $pwd/interproscan/$input[$#input]-interproscanRC2 -dp -goterms -iprlookup -pa -f TSV,XML,GFF3\""); # v1.01 udpate
		# For the above one to be used I should fix phobius.pl for phobius -- done Fixed!
		# ncoils for coils -- done Fixed!
		# blastall for PIRSF, Panther, ProDom -- done Fixed!
		
		# So modified to test the installation 
		#push(@commands,"bsub -q $queue -K -J iprscan5_$count -o $pwd/interproscan/out_interproscan_rc2_$input[$#input]-$count -R \"rusage[mem=4096]\" \"source interproscan-5;interproscan.sh -appl TIGRFAM,PrositeProfiles,PrositePatterns,PRINTS,SuperFamily,Gene3d,PfamA,TMHMM,HAMAP -i $_ -b $pwd/interproscan/$input[$#input]-interproscanRC2 -dp -goterms -iprlookup -pa -f TSV\"");
		# for testing only
		#push(@commands,"bsub -q $queue -K -J iprscan5_$count -o $pwd/interproscan/out_interproscan_rc2_$input[$#input]-$count -R \"rusage[mem=4096]\" \"source interproscan-5;interproscan.sh -appl Coils -i $_ -b $pwd/interproscan/$input[$#input]-interproscanRC2 -dp -goterms -iprlookup -pa -f TSV\" && echo INTERPROSCAN for chunk $_ done");
		$count++;
	}
	&process_cmd(\@commands,"INTERPROSCAN");

	# Now concatenate the results 
	# add line to end of the file (http://www.thegeekstuff.com/2009/11/unix-sed-tutorial-append-insert-replace-and-count-file-lines/)
	@commands=(); # empty the interproscan commands
	# v1.01 update
	push(@commands,"bsub -q $queue -K -J cat_interpro_tsv -o $pwd/interproscan/out_cat_interproscan_tsv \"cat $pwd/interproscan/$input*-interproscanRC2.tsv | sort -k1,1 | sed '\$ a END' > $pwd/combined_$input.iprscan.tsv\"");
	push(@commands,"bsub -q $queue -K -J cat_interpro_xml -o $pwd/interproscan/out_cat_interproscan_xml \"cat $pwd/interproscan/$input*-interproscanRC2.xml > $pwd/combined_$input.iprscan.xml\"");
	push(@commands,"bsub -q $queue -K -J cat_interpro_gff3 -o $pwd/interproscan/out_cat_interproscan_gff3 \"cat $pwd/interproscan/$input*-interproscanRC2.gff3 > $pwd/combined_$input.iprscan.gff3\"");

	&process_cmd(\@commands,"INTERPROSCAN");
	#print "Interproscan jobs executed!!\n";
	return 1;
} #end of subroutine

#########################
###Blast2GO##############
#########################

sub blast2go() {
	chomp(my @files=`ls $pwd/leaff/*.txt`);
	my $count = 1;
	my @commands=();
	foreach (@files) { # if given nucleotide change blastp to blastx
		my @input=split (/\//);
		
		# /tgac/references/databases/blast/nr -- new
		# /data/references/databases/blast/nr -- old
		
		# For input peptides only "blastp"
		push(@commands,"bsub -q $queue -K -J Blast2GOblastp_$count -o $pwd/blastpForBlast2go/out_blastpForBlast2GO_$input[$#input]-$count -n 4 -R \"span[ptile=4] rusage[mem=5120]\" \"source jre-7.11;source blast-2.2.22;blastall -p blastp -a 4 -d /tgac/public/databases/blast/ncbi/nr -i $_ -e 1e-4 -b 50 -v 50 -m 7 -o $pwd/blastpForBlast2go/$input[$#input]-vs-nr_blastpForBlast2GO.xml\"");	# v1.01 update to blast db location
		$count++;
	}
	&process_cmd(\@commands,"BlastpForBlast2GO");
	#print "Blastp jobs for Blast2go executed!!\n";

	# Now concatenate the results
	@commands=(); # empty the blastp commands
	push(@commands,"bsub -q $queue -K -J cat_Blast2GOblastp -o $pwd/blastpForBlast2go/out_cat_blastpForBlast2GO \"cat $pwd/blastpForBlast2go/*blastpForBlast2GO.xml > $pwd/combined_$input.blastpForBlast2GO.xml\"");
	&process_cmd(\@commands,"BlastpForBlast2GO");

	@commands=(); # empty the cat commands
	if(!-e "$pwd/combined_$input.blastpForBlast2GO.xml"){
		die "Exiting .. blastpForBlast2GO didn't generate any output for Blast2GO\n";
		exit;
	}
	else {
		push(@commands,"bsub -q $queue -K -o $pwd/blast2go/out_blast2go -J blast2go -n 8 -R \"span[ptile=8] rusage[mem=10240]\" \"source x11-7.0.0.22;source blast2go-2.5.0;java -cp /tgac/software/testing/blast2go/2.5.0/src/b2g4pipe/*:/tgac/software/testing/blast2go/2.5.0/src/b2g4pipe/ext/*: es.blast2go.prog.B2GAnnotPipe -in $pwd/combined_$input.blastpForBlast2GO.xml -out $pwd/blast2go/combined_$input.blastpForBlast2GO.xml.blast2go -prop /tgac/software/testing/blast2go/2.5.0/src/b2g4pipe/b2gPipe.properties -v -annot -annex -goslim -dat -img\"");
		&process_cmd(\@commands,"BLAST2GO");
	}
	#print "Blast2go job executed!!\n";
	return 1;
} #end of subroutine

#########################
###PARSE#################
#########################
sub parse {

	my @commands=();

	# Clean up first
	push(@commands, "cat $pwd/leaff/out_leaff_* > $pwd/leaff/out_leaff;rm $pwd/leaff/out_leaff_*;mv $pwd/error.log $pwd/blastpForBlast2go/blastp_error.log;cat $pwd/blastpForBlast2go/out_blastpForBlast2GO_* > $pwd/blastpForBlast2go/out_blastpForBlast2GO;rm $pwd/blastpForBlast2go/out_blastpForBlast2GO_*;cat $pwd/interproscan/out_interproscan_rc2_* > $pwd/interproscan/out_interproscan; rm $pwd/interproscan/out_interproscan_rc2_*;mv $pwd/temp $pwd/interproscan_temp ");
	&process_cmd(\@commands,"Cleanup");

	#interproscan	
	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"parse_interproscan_tsv.pl $pwd/combined_$input.iprscan.tsv > $pwd/temp_folder/combined_$input.iprscan.tsv.out\"");
	&process_cmd(\@commands,"INTERPROSCAN");

	#blast2go
	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"parse_blast2go_annot.pl $pwd/blast2go/combined_$input.blastpForBlast2GO.xml.blast2go.annot > $pwd/temp_folder/combined_$input.blastpForBlast2GO.xml.blast2go.annot.out\"");
	&process_cmd(\@commands,"BLAST2GO");
	
	#combine results
	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"cat $pwd/temp_folder/combined_$input.blastpForBlast2GO.xml.blast2go.annot.out $pwd/temp_folder/combined_$input.iprscan.tsv.out | sort -k1,1 > $pwd/temp_folder/combined_Blast2GO_interproscan.results.txt\"");
	&process_cmd(\@commands,"Parse");
	
	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"parse_combined_interproscan_blast2go.pl $pwd/temp_folder/combined_Blast2GO_interproscan.results.txt > $pwd/temp_folder/combined_Blast2GO_interproscan.results.txt.out\"");
	&process_cmd(\@commands,"Parse");

	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"grep '^>' $pwd/$input | sed 's/^>//g' | awk '{print \\\$1\\\"\\t\\\"\\\$1}' > $pwd/temp_folder/$input-contigID.txt\"");
	&process_cmd(\@commands,"Parse");

	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"parse_compareIDs_input_interproscan_blast2go.pl $pwd/temp_folder/$input-contigID.txt $pwd/temp_folder/combined_Blast2GO_interproscan.results.txt.out | sort -k1,1 > $pwd/$input-annotation.tsv\"");
	&process_cmd(\@commands,"Parse");

	@commands=();
	push(@commands, "bsub -q $queue -K -o $pwd/temp_folder/out_parse -J parse \"sed -i \'1i #ID\\tBlast2GO_GO_term\\tBlast2GO_EC_Number\\tBlast2GO_GO_Description\\tInterproscan_IPR\\tInterproscan_IPR_Description\\tInterproscan_GO_term\\tInterproscan_EC_Number\\tInterproscan_Pathways\' $pwd/$input-annotation.tsv\"");
	&process_cmd(\@commands,"Parse");
}

#########################
###RUN COMMANDS##########
#########################

sub process_cmd {
	my ($commands,$msg) = @_;
	
	if ($msg) {
        print "\n";
        print "#######################################\n";
        print "##    Starting ... $msg\n";
        print "#######################################\n";
    }
	my $start_time = time();
	my @threads;
	foreach my $cmd (@$commands) {
		#my $t = threads->new('system_call',$cmd);
		my $t = threads->new('system_call',$cmd,$msg);
		push(@threads,$t);
	}
	foreach (@threads) {
		$_->join;
	}
	my $end_time = time();
	print "CMD finished : $msg (" . ($end_time - $start_time) . " seconds)\n";
}

#	get ideas from here : http://stackoverflow.com/questions/1380516/how-can-i-use-perls-system-call-to-spawn-independent-threads
#	and here
#	https://wiki.bc.net/atl-conf/pages/viewpage.action?pageId=20548191

sub system_call {
#	my $cmd=shift;
#	system("date");
#	print "$cmd\n";
#	system($cmd);           # || use @_ instead of $cmd if not collecting to variable

	my @cmd=@_;     # $cmd[0] is the command $cmd[1] is the message
	#print "$cmd[0]------1\n$cmd[1]------2\n";
	system("date >> $cmd[1].log 2>&1");
	open (MYFILE, '>>', "$cmd[1].log");
	print MYFILE "$cmd[0]\n";
	system("$cmd[0] >> $cmd[1].log 2>&1");           # || use @_ instead of $cmd if not collecting to variable

}

#########################
###MAIN##################
#########################
sub main() {
#	http://search.cpan.org/~szabgab/Parallel-ForkManager-1.03/lib/Parallel/ForkManager.pm
	Parallel::ForkManager->import(); ## added new 06/12/2013
	my $pm = Parallel::ForkManager->new($THREADS_NUM);
	{
		$pm->start and next; # do the fork
		&leaff();				# Run leaff
		$pm->finish; # do the exit in the child process
	}
	$pm->wait_all_children;
	#print "Done! Splitting fasta\n";

#	http://perldoc.perl.org/threads.html
	my $thr1 = threads->create('interproscanRC2');
	my $thr2 = threads->create('blast2go');
	$thr1->join;
	$thr2->join;
	#print("Done! Interproscan and Blast2GO_blast\n");
	
	my $pm2 = Parallel::ForkManager->new($THREADS_NUM);
	{
		$pm2->start and next; # do the fork
		&parse();
		$pm2->finish; # do the exit in the child process
	}
	$pm2->wait_all_children;
	#print "Done! Parsing\n";

	print "\n\n";
	print "#################################################################\n";
    print "The Functional annotation files are below:\n$pwd/$input-annotation.tsv\n";
    print "\nRaw Interproscan result can be found here:\n";
    print " -	TSV format : $pwd/combined_$input.iprscan.tsv\n";
    print " -	XML format : $pwd/combined_$input.iprscan.xml\n";
    print " -	GFF3 format : $pwd/combined_$input.iprscan.gff3\n";
    print "\nRaw Blast xml ouput file used for Blast2GO can be found here:\n$pwd/combined_$input.blastpForBlast2GO.xml\n";
    print "\nRaw Blast2GO results can be found here:\n$pwd/blast2go/combined_$input.blastpForBlast2GO.xml.blast2go*\n";
    print "#################################################################\n";
	exit;
}
main();
# end

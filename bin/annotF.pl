#!/usr/bin/env perl

#	This pipeline with execute interproscan, blastp for blast2go for input "peptides"
#	at the same time and finally combine all the results as a functional annotation file

#	This script use packages : LEAFF, INTERPRO, BLAST2GO (including blastp)

# 	Version 1.02
#	- Adding option to run on SLURM environment
#	- Added option to use blastdb from user # Monday, 12 June 2017, 12:05PM
# 	- # Tuesday, 22 January 2019, 12:07PM
# 	- Increased interproscan memory from 8GB to 20GB
# 	- Changed occurences of tgac to ei


# 	Version 1.01
#	- Added XML,GFF3 interproscan output formats
#	- Changed the blast database location from this /ei/references/databases/blast/nr to this /ei/public/databases/blast/ncbi/nr

# 	Version 1.0
#	- First release

# AUTHOR: Gemy George Kaithakottil (Gemy.Kaithakottil@earlham.ac.uk || Gemy.Kaithakottil@gmail.com)


use strict;
use warnings;
use File::Basename;
use File::Path qw(make_path remove_tree);
use Parallel::ForkManager;
#use Getopt::Long qw(:config no_ignore_case pass_through); 
use Getopt::Long; 
use Cwd;
use FindBin;
use lib "$FindBin::RealBin/PerlLib";
use HPC::GridRunner;

my $help_flag;
my $input;
my $nr; # Monday, 12 June 2017, 12:05PM
my $jobScheduler;
my $num = 1000;
my $queue = "ei-medium";
my $prop = "/ei/software/testing/blast2go/2.5.0/src/b2g4pipe/b2gPipe.30Jan2017_database.properties";
my $help;
my $prog = basename($0);

my $usage = <<_EOUSAGE_;

#########################################################################################################
#
#    AnnotF (v1.02) is a functional annotation pipeline for the annotation of query proteins
#   
#        Usage: $prog [options] --fasta <[protein.fa || /path/to/protein.fa]>
#
#	Required:
#	--fasta <string>            fasta file, properly formatted 
#	                            (*IMPORTANT*:Please read fasta file requirements below)
#	--nr <string>               provide path to blast nr database
#	                            /path/to/blast/nr
#	                            Eg: /ei/public/databases/blast/ncbi/nr_20170116/nr
#
#	General Options:
#	--chunk_size <int>          required chunk size (fasta per chunk) for faster execution of jobs
#	                            DO NOT exceed 1000 fasta sequences per chunk or fasta file
#	                            (Default: $num, fasta sequences per chunked fasta file)
#	--queue <string>            queue to submit jobs
#	                            (Default: ei-medium, for SLURM)
#	--prop <property-file>      file containing the Blast2GO application settings 
#	                            (Default: /ei/software/testing/blast2go/2.5.0/src/b2g4pipe/b2gPipe.30Jan2017_database.properties)
#	--help                      print this option menu and quit
#							
#    NOTE: 
#    *Fasta file requirements:
#    - There should not be any "dots" in sequences or "empty/blank lines" in the fasta file
#    - Only fasta header should be present in fasta sequence and no pipe ("|") character is
#      allowed in fasta header name as it will cause issues with InterproScan. 
#      For example, a fasta sequence name header ">gb|xyz|abc" is NOT allowed, since it has pipe character.
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
#    Example SLURM command:
#    sbatch -p ei-medium -J annotF -o out_annotF.%j.log -c 1 --mem 5G --wrap "source annotF-1.02 && annotF --queue ei-medium --fasta test_proteins.fasta --nr /ei/public/databases/blast/ncbi/nr_20170116/nr"
#
#
# Contact : Gemy George Kaithakottil (Gemy.Kaithakottil\@earlham.ac.uk || Gemy.Kaithakottil\@gmail.com)
#########################################################################################################
_EOUSAGE_
    ;

unless (@ARGV) {
    die "$usage\n";
}

&GetOptions ( 'fasta=s' => \$input,
			  # 'jobScheduler=s' => \$jobScheduler,
			  'nr=s' => \$nr,
              'chunk_size=i' => \$num,
              'queue=s' => \$queue,
              'prop=s' => \$prop,
              'h|help' => \$help,

);

#if ($help_flag || !$input || !$num) {
if ($help) {
    die $usage;
}

if (@ARGV) {
    die "# ERROR, do not understand options: @ARGV\n";
}

unless ($input) {
    die "# ERROR: No input fasta provided.\n$usage\n";
}

unless(-e $input) {
	die "# ERROR: Cannot open the input file - $input\n$!\n";
}

# Link the protein file to the current location
my $input_basename = basename($input);
qx(ln -s $input .) unless (-e $input_basename);
$input = $input_basename;


unless ($nr) {
    die "# ERROR: No path to nr database provided.\n$usage\n";
}

# This is going to be the file prefix used at any stage of the pipeline
my $prefix = "$input.chunk";

# Get current working directory
my $pwd = getcwd;

#To pass number of chunks created in leaff step;
my $chunk=0;

# CONSTANTS
# Number of jobs to be executed at a time -- default 1000 jobs
my $THREADS_NUM = 1000; # http://bickson.blogspot.co.uk/2011/12/multicore-parser-part-2-parallel-perl.html

# Print version status:
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
print "%%                            annotF - v1.02                            %%\n";
print "%%   Gemy.Kaithakottil\@earlham.ac.uk || Gemy.Kaithakottil\@gmail.com     %%\n";
print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";

# Create directories
print "\n";
print "#######################################\n";
print "##    Creating necessary directories..\n";
print "#######################################\n";
my $leaff_dir = "$pwd/1.leaff";
my $interproscan_dir = "$pwd/2.interproscan";
my $blastp_dir = "$pwd/3.blastp";
my $blast2go_dir = "$pwd/4.blast2go";
my $temp_folder_dir = "$pwd/5.temp_folder";
make_path($leaff_dir, $interproscan_dir, $blastp_dir, $blast2go_dir, $temp_folder_dir, {verbose => 1});


# DEFAULT BLAST DB LOCATION
# my $blast_DB="/ei/public/databases/blast/ncbi/nr";
my $blast_DB = $nr; # Monday, 12 June 2017, 12:04PM

# create configuration file as required
sub write_conf {
	my $filename = shift;
	my $queue = shift;
	my $ncpus = shift;
	my $mem = shift;
	my $cmds_per_node = shift;
	open(my $fh, '>', $filename) or die "Could not write file '$filename' $!\n";
	print $fh "# grid type\n";
	print $fh "grid=SLURM\n";
	print $fh "# template for a grid submission\n";
	print $fh "cmd=sbatch -p $queue -c $ncpus --mem=$mem\n";
	print $fh "# number of grid submissions to be maintained at steady state by the submission system\n";
	print $fh "max_nodes=1000\n";
	print $fh "# number of commands that are batched into a single grid submission job.\n";
	print $fh "cmds_per_node=$cmds_per_node\n";  # need to change it later
	close $fh;
}

# write output files locally
# 	- input : location with prefix, with wild characters
# 	- output : output filename with absolute path
sub write_files {
	my $input_files = shift;
	my $output_filename = shift;

	chomp(my @files= qx(ls $input_files));
	open OUTPUT_FILE, ">", $output_filename or die "# Fatal error: Cannot write file '$output_filename'. $!\n";
	foreach my $each_file (@files) {
		if (open my $in, '<', $each_file) {
			while (my $line = <$in>) {
		    	print OUTPUT_FILE $line;
		}
		close $in;
		} else {
			warn "Could not open '$each_file' for reading\n";
		}
	}
	close(OUTPUT_FILE);
}


#########################
##LEAFF##################
#########################
sub leaff() {
	#print $log "Inside leaff\n";
	my $file = shift;
	my $conf = "$leaff_dir/SLURM.conf";
	my $ncpus = 1;
	my $mem = 2048;
	my $cmds_per_node = 1;
	write_conf($conf,$queue,$ncpus,$mem,$cmds_per_node);

	my $count=0;
	my $lines;
	$lines = qx(grep -c '^>' $file); # MAIN
	chomp($lines);
	my $total = ($lines/$num);
	my @commands=();
	$num = $num-1;
	for (my $i=0;$i<=($lines-1);$i++) {
		$count++;
		my $j=$i+$num; # like 9,19,29...
		if ($j <= $lines && ($j != $lines)) {
			#print $log "Starting ",($i+1),"-",($j+1)," End\n";
			push (@commands, "leaff -f $file -S $i $j > $leaff_dir/$prefix-$count.txt");
			$i=$j;
		}
		else {
			$j = ($lines-1);
			#print $log "Starting ",($i+1),"-",($j+1)," End**Last File**\n";
			push (@commands, "leaff -f $file -S $i $j > $leaff_dir/$prefix-$count.txt");
			last;
		}
	}
	$chunk = $count; # see whether its used again or else delete it
	
	#http://stackoverflow.com/questions/3991143/how-can-i-pass-two-arrays-and-a-string-to-a-perl-subroutine
	#pass array as reference
	&process_cmd(\@commands,"LEAFF",$conf);
} #end of subroutine

#########################
###Interproscan_rc2######
#########################

sub interproscan() {

	# write configuration file 
	my $file = $input;
	my $conf = "$interproscan_dir/SLURM.conf";
	my $ncpus = 6;
	# my $mem = 8192;
	my $mem = 20480;
	my $cmds_per_node = 1;
	write_conf($conf,$queue,$ncpus,$mem,$cmds_per_node);

	chomp(my @files= qx(ls $leaff_dir/$prefix*.txt)); # read from leaff directory

	my $count = 1;
	my @commands=();
	foreach my $filepath (@files) {
		my $filebase = basename ($filepath);
		push(@commands,"interproscan.sh -i $filepath -b $interproscan_dir/$filebase-iprscan5 -dp -goterms -iprlookup -pa -f TSV,XML,GFF3"); # v1.02 udpate
	}
	&process_cmd(\@commands,"INTERPROSCAN",$conf);

	# RUN LOCALLY
	print "#INFO: Concatenate interproscan files\n";

	# write TSV files
	write_files("$interproscan_dir/$prefix*.tsv","$pwd/$input.interproscan.tsv");

	system("cp -a $pwd/$input.interproscan.tsv $temp_folder_dir/$input.interproscan.tsv.bkp && cat $temp_folder_dir/$input.interproscan.tsv.bkp | sort -k1,1 | sed '\$ a END' > $temp_folder_dir/$input.interproscan.tsv; rm $temp_folder_dir/$input.interproscan.tsv.bkp");

	# write XML files
	write_files("$interproscan_dir/$prefix*.xml","$pwd/$input.interproscan.xml");

	# write GFF3 files
	write_files("$interproscan_dir/$prefix*.gff3","$pwd/$input.interproscan.gff3");

} #end of subroutine

#########################
###Blastp################
#########################

sub blastp() {

	# write configuration file 
	my $file = $input;
	my $conf = "$blastp_dir/SLURM.conf";
	my $ncpus = 4;
	my $mem = 40960;
	my $cmds_per_node = 1;
	write_conf($conf,$queue,$ncpus,$mem,$cmds_per_node);

	chomp(my @files=qx(ls $leaff_dir/$prefix*.txt)); # read from leaff directory);
	my $count = 1;
	my @commands=();

	foreach my $filepath (@files) {
		my $filebase = basename ($filepath);
		# For input peptides only "blastp"
		push(@commands,"blastall -p blastp -a $ncpus -d $blast_DB -i $filepath -e 1e-4 -b 50 -v 50 -m 7 -o $blastp_dir/$filebase.blastp.xml");
		$count++;
	}
	&process_cmd(\@commands,"BLASTP",$conf);
	#print "Blastp jobs for Blast2go executed!!\n";

	# Now concatenate the results
	@commands=(); # empty the blastp commands

	# RUN LOCALLY
	print "#INFO: Concatenate blastp files\n";

	# write TSV files
	write_files("$blastp_dir/$prefix*.blastp.xml","$pwd/$input.blastp.xml");
} #end of subroutine



#########################
###Blast2GO##############
#########################

sub blast2go() {

	# make sure that the blastp output file is present
	if(!-e "$pwd/$input.blastp.xml"){
		die "Exiting .. blastp results for Blast2GO does not exist. Try rerunning annotF..\n";
		exit;
	}
	else {
		# write configuration file 
		my $file = $input;
		my $conf = "$blast2go_dir/SLURM.conf";
		my $ncpus = 8;
		my $mem = 20480;
		my $cmds_per_node = 1;
		write_conf($conf,$queue,$ncpus,$mem,$cmds_per_node);
		my @commands=();
		my $java_Xmx=$mem . "m";

		die "# Fatal error: Cannot file $prop file. Please provide correct link to properties file and re-run\n" unless (-e $prop);

		push(@commands,"source blast2go-2.5.0; unset DISPLAY; cp -a $prop $blast2go_dir/b2gPipe.properties; java -Xmx$java_Xmx -cp /ei/software/testing/blast2go/2.5.0/src/b2g4pipe/*:/ei/software/testing/blast2go/2.5.0/src/b2g4pipe/ext/*: es.blast2go.prog.B2GAnnotPipe -in $pwd/$input.blastp.xml -out $blast2go_dir/$input.blast2go -prop $blast2go_dir/b2gPipe.properties -v -annot -annex -goslim -dat -img");
		&process_cmd(\@commands,"BLAST2GO",$conf);
	}
} #end of subroutine


#########################
###PARSE#################
#########################
sub parse {

	my @commands=();

	# write configuration file 
	my $file = $input;
	my $conf = "$temp_folder_dir/SLURM.conf";
	my $ncpus = 2;
	my $mem = 5120;
	my $cmds_per_node = 1;
	write_conf($conf,$queue,$ncpus,$mem,$cmds_per_node);

	#interproscan	
	@commands=();
	push(@commands, "parse_interproscan_tsv.pl $temp_folder_dir/$input.interproscan.tsv > $temp_folder_dir/$input.interproscan.tsv.out && parse_blast2go_annot.pl $blast2go_dir/$input.blast2go.annot > $temp_folder_dir/$input.blast2go.annot.out && cat $temp_folder_dir/$input.blast2go.annot.out $temp_folder_dir/$input.interproscan.tsv.out | sort -k1,1 > $temp_folder_dir/$input.blast2go_interproscan.results.txt && parse_combined_interproscan_blast2go.pl $temp_folder_dir/$input.blast2go_interproscan.results.txt > $temp_folder_dir/$input.blast2go_interproscan.results.txt.out && awk '/^>/ {split(\$0,a,\" \");print a[1]\"\\t\"a[1];}' $input | sed 's/>//g' > $temp_folder_dir/$input-contigID.txt && parse_compareIDs_input_interproscan_blast2go.pl $temp_folder_dir/$input-contigID.txt $temp_folder_dir/$input.blast2go_interproscan.results.txt.out | sort -k1,1V > $pwd/$input.annotF-annotation.tsv && sed -i \'1i #ID\\tBlast2GO_GO_term\\tBlast2GO_EC_Number\\tBlast2GO_GO_Description\\tInterproscan_IPR\\tInterproscan_IPR_Description\\tInterproscan_GO_term\\tInterproscan_EC_Number\\tInterproscan_Pathways\' $pwd/$input.annotF-annotation.tsv");
	&process_cmd(\@commands,"Generate_output_files",$conf);

}

#########################
###RUN COMMANDS##########
#########################

sub process_cmd {
	my ($commands,$msg,$conf) = @_;
	
	if ($msg) {
        print "\n";
        print "#######################################\n";
        print "##    Starting ... $msg\n";
        print "#######################################\n";
    }
	my $start_time = time();
	my $grid_runner = new HPC::GridRunner($conf, "$msg.htc_cache_success");
    my $ret = $grid_runner->run_on_grid(@$commands);
    if ($ret) {
        die "Error, not all $msg commands completed successfully.  Cannot continue.\n";
    }
	my $end_time = time();
	print "CMD finished : $msg (" . ($end_time - $start_time) . " seconds)\n";
}

#	get ideas from here : http://stackoverflow.com/questions/1380516/how-can-i-use-perls-system-call-to-spawn-independent-threads
#	and here
#	https://wiki.bc.net/atl-conf/pages/viewpage.action?pageId=20548191

#########################
###MAIN##################
#########################
sub main() {

	# Check that the environment have the executable variables
	my @tools = qw (leaff interproscan.sh java blastall parse_interproscan_tsv parse_blast2go_annot parse_combined_interproscan_blast2go parse_compareIDs_input_interproscan_blast2go);
	foreach my $tool (@tools) {
		my $tool_path = qx(which $tool);
		chomp ($tool_path);
		die "\n## Fatal error: No '$tool' command available in the env PATH. Please make sure you have the executables before executing the pipeline\n" unless ( $tool_path );
	}

	# http://search.cpan.org/~szabgab/Parallel-ForkManager-1.03/lib/Parallel/ForkManager.pm
	Parallel::ForkManager->import(); ## added new 06/12/2013
	my $pm = Parallel::ForkManager->new($THREADS_NUM);
	{
		$pm->start and next; # do the fork
		&leaff($input);				# Run leaff
		$pm->finish; # do the exit in the child process
	}
	$pm->wait_all_children;
	#print "Done! Splitting fasta\n";

	my $main_pm = Parallel::ForkManager->new($THREADS_NUM);
	$main_pm->run_on_start(
		sub { my ($pid,$ident)=@_;
		# print "** $ident started, pid: $pid\n";
		&interproscan;
		}
	);

	# get it's exit code
	$main_pm->run_on_finish(
		sub { my ($pid, $exit_code, $ident) = @_;
		# print "** $ident just got out of the pool with PID $pid and exit code: $exit_code\n";
		&blast2go;
		}
	);

	{
		my $pid = $main_pm->start("BLASTP") and next; # do the fork
		&blastp;
		$main_pm->finish; # do the exit in the child process
	}
	print STDERR "Waiting for Children...\n";
	$main_pm->wait_all_children;
	print STDERR "Everybody is out of the pool!\n";


	my $pm2 = Parallel::ForkManager->new($THREADS_NUM);
	{
		$pm2->start and next; # do the fork
		&parse();
		$pm2->finish; # do the exit in the child process
	}
	$pm2->wait_all_children;
	#print "Done! Parsing\n";


	# Clean directories
	print "\n";
	print "#######################################\n";
	print "##    Clean the output directory\n";
	print "#######################################\n";
	my @list=();
	# http://stackoverflow.com/questions/11653762/check-for-existence-of-directory-in-perl-with-wildcard
	@list = `ls -d $pwd/farmit.*` if (grep -d, glob "$pwd/farmit.*");
	push (@list, `ls -d $pwd/temp`) if (-d "$pwd/temp");
	chomp @list;
	if (@list){
		foreach my $dir (@list) {
			print $dir,"\n";
			remove_tree($dir);
		}
	}


	print "\n\n";
	print "#################################################################\n";
	print "#################################################################\n";
	print "#################################################################\n";
    print "The functional annotation files are below:\n";
    print "\nAnnotF output:\n";
    print " -	TSV format : $pwd/$input.annotF-annotation.tsv\n";
    print "\nRaw InterProScan v5.22-61.0 result can be found here:\n";
    print " -	TSV format : $pwd/$input.interproscan.tsv\n";
    print " -	XML format : $pwd/$input.interproscan.xml\n";
    print " -	GFF3 format : $pwd/$input.interproscan.gff3\n";
    print "\nRaw Blast xml ouput file used for Blast2GO can be found here:\n";
    print " -	XML format : $pwd/$input.blastp.xml\n";
    print "\nRaw Blast2GO v2.5.0 (Database 30Jan2017) results can be found here:\n";
    print " -	$blast2go_dir/$input.blast2go*\n";
    print "#################################################################\n";
    print "#################################################################\n";
    print "#################################################################\n";



	exit;
}
main();
# end

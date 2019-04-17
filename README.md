# AnnotF - A functional Annotation pipeline

AnnotF is a pipeline for the functional annotation of peptide sequences while making use 
of commonly used functional annotation tools - Interproscan, Blast and Blast2GO. The pipeline 
executes jobs in parallel and it provide users a consolidated tab-delimited output along with 
individual tools output.

## Directed acyclic graph (DAG)
The below graph explains the basic workflow 

## Requirements
```
perl > 5.16.0
  Parallel::ForkManager
LEAFF - http://kmer.sourceforge.net/wiki/index.php/LEAFF_User%27s_Guide
blast-2.2.22 - ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.22/
blast2go-2.5.0 - https://www.blast2go.com/support/blog/22-blast2goblog/110-local-blast2go-database-installation
interproscan-5.22.61 - https://www.ebi.ac.uk/about/news/service-news/InterPro-61.0
```

## Run AnnotF
AnnotF requires only one mandatory argument. The query peptide file, for example a protein fasta file is below:

CAUTION: 
Only fasta header should be present and no pipe ("|") allowed in fasta header as it will cause issues with InterproScan.
```
>A2YIW7
MAAEEGVVIACHNKDEFDAQMTKAKEAGKVVIIDFTASWCGPCRFIAPVFAEYAKKFPGAVFLKVDVDELKEVAEKYNVE
AMPTFLFIKDGAEADKVVGARKDDLQNTIVKHVGATAASASA
>comp21620_c0_seq1
MMGKGKVLVYHFLPSIFLVNITYYIYIYIYIYIFLYIHTYNYIHIYILLMASITGSAVSI
SSFSCSFKLNQASARVSTLNSVPFSINGKSFPSIKLRPAQRFQVSCMAAKPETVEKVCGI
VRKQLAIAADTEITGESKFAALGADSLDTVEIVMGLEEEFGISVEEESAQTIATVQDAAD
LIEKLLA
```
Example commands on SLURM:
```
sbatch -p ei-medium -J annotF -o out_annotF.%j.out -e out_annotF.%j.err -c 1 --mem 5G --wrap "source annotF-1.02; annotF --queue ei-medium --fasta test_proteins.fasta --nr /tgac/public/databases/blast/ncbi/nr_20170116/nr"
```

INFO:
If your jobs have been killed due to HPC failure/issues, just re-lauch the intial command you have used and annotF should pick up from where it stopped.


You cound find the help menu by running "annotF -h" in the terminal.
```
#########################################################################################################
#
#    AnnotF (v1.02) is a functional annotation pipeline for the annotation of query proteins
#
#        Usage: annotF [options] --fasta <[protein.fa || /path/to/protein.fa]>
#
#	Required:
#	--fasta <string>            fasta file, properly formatted
#	                            (*IMPORTANT*:Please read fasta file requirements below)
#	--nr <string>               provide path to blast nr database
#	                            /path/to/blast/nr
#	                            Eg: /tgac/public/databases/blast/ncbi/nr_20170116/nr
#
#	General Options:
#	--chunk_size <int>          required chunk size (fasta per chunk) for faster execution of jobs
#	                            DO NOT exceed 1000 fasta sequences per chunk or fasta file
#	                            (Default: 1000, fasta sequences per chunked fasta file)
#	--queue <string>            queue to submit jobs
#	                            (Default: tgac-medium, for SLURM)
#	--prop <property-file>      file containing the Blast2GO application settings
#	                            (Default: /tgac/software/testing/blast2go/2.5.0/src/b2g4pipe/b2gPipe.30Jan2017_database.properties)
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
#    sbatch -p tgac-medium -J annotF -o out_annotF.%j.out -e out_annotF.%j.err -c 1 --mem 5G --wrap "source annotF-1.02; annotF --queue tgac-medium --fasta test_proteins.fasta --nr /tgac/public/databases/blast/ncbi/nr_20170116/nr"
#
#
# Contact : Gemy George Kaithakottil (Gemy.Kaithakottil@earlham.ac.uk || Gemy.Kaithakottil@gmail.com)
#########################################################################################################
```
## Example

Test file is provided under annotF home bin folder:

Sample Input file:
```
examples/test_proteins.fasta
```
You could copy the file across to your working directory and run annotF and compare your output to the provided example output to check whether all packages work as it should.

Sample output file is provided here:
```
examples/test_proteins.fasta.annotF-annotation.tsv
```

## Results

Complete list of output files present in the output directory includes:

AnnotF output
-	TSV format : <input_file>.annotF-annotation.tsv

Raw InterProScan v5.22-61.0 result can be found here:
-	TSV format : <input_file>.interproscan.tsv
-	XML format : <input_file>.interproscan.xml
-	GFF3 format : <input_file>.interproscan.gff3

Raw Blast xml ouput file used for Blast2GO can be found here:
-	XML format : <input_file>.blastp.xml

Raw Blast2GO v2.5.0 (Database 30Jan2017) results can be found here:
-	4.blast2go/<input_file>.blast2go*


## Contact

If you encounter a bug or if you have any suggestions please contact me through email below with subject "annotF-<query>".

Gemy.Kaithakottil@earlham.ac.uk
Gemy.Kaithakottil@gmail.com

## Acknowledgements

I would like to thank my group leader David.Swarbreck@earlham.ac.uk for his support and giving me the opportunity to develop AnnotF.

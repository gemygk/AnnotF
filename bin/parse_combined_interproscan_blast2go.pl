#!/usr/bin/env perl


#	parse the iprscan output file

# 	Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk)
# 	Version 0.1
#	Part of interproscan and blast2go pipeline

use strict;
#use warnings;

my $usage = "\nUsage: perl $0 <final.txt>

This input final.txt file is a combination of output from running the scripts blast2goAnnocationParse_pipeline.pl (on blast2go annot file) and iprscan_parse_rc2_compressed.pl (on iprscan_rc2 .tsv file)
Note: The file needed to be sorted by first column (sort -k1,1 <filename>, will do)

";
my $file = shift or die $usage;
my $inf;
my $true=1;
my $for_pass="";
my ($id,$go_ipr,$ec_iprid,$godesc_iprdesc,$iprdesc,$gotag,$iprgo,$pathtag,$pathways);
my ($nid,$ngo_ipr,$nec_iprid,$ngodesc_iprdesc,$niprdesc,$ngotag,$niprgo,$npathtag,$npathways);
my (@blast2go_go,@iprscan_ipr,@iprscan_go,@iprscan_iprdesc,@iprscan_path)=();
my (%blast2go_goH,%iprscan_iprH,%iprscan_goH,%iprscan_iprdescH,%iprscan_pathH,%iprscan_ecH)=();
my (@blast2go_goA,@iprscan_iprA,@iprscan_goA,@iprscan_iprdescA,@iprscan_pathA,@iprscan_ecA)=();
my $count=0;
my $rare = 1;
my $boolean = 0;
my $dummy;
my $allright=1;
my $notallright=0;
my $newline;
my $gone=0;
open $inf, "<" . $file;

#print the header
#print "ID\tBlast2GO_GO_term\tBlast2GO_EC_Number\tBlast2GO_GO_Description\tInterproscan_IPR\tInterproscan_IPR_Description\tInterproscan_GO_term\tInterproscan_EC_Number\tInterproscan_Pathways\n";
while (<$inf>) {
	chomp();
	if($boolean) {
		$dummy = $_;
                $_ = $newline;
                $newline = $dummy;
                $boolean = 0;
                $rare = 0;
	}
	($id,$go_ipr,$ec_iprid,$godesc_iprdesc,$iprdesc,$gotag,$iprgo,$pathtag,$pathways) = split(/\t/);
	if($rare) { #true first time
    	$newline = <$inf>;
    	chomp($newline);
    }
	($nid,$ngo_ipr,$nec_iprid,$ngodesc_iprdesc,$niprdesc,$ngotag,$niprgo,$npathtag,$npathways) = split(/\t/, $newline);
	$rare = 1;
	if ($id eq $nid) {
		
		# identify blast2go GO repeats and remove
		if ($go_ipr =~ /^GO/) {
                       if ($go_ipr =~ m/\|/) {
                               @blast2go_go = split(/\|/,$go_ipr);
                       }
                       else {
                               push(@blast2go_go,$go_ipr);
                       }
		}
		foreach (@blast2go_go) {
			unless (exists ($blast2go_goH{$_})) {
				push (@blast2go_goA,$_);
				$blast2go_goH{$_} += 1;
			}
		}
		
		# identify iprscan IPR id repeats and remove
		if ($ngo_ipr =~ /^#IPRID#/) {
                       if ($nec_iprid =~ m/\|/) {
                               @iprscan_ipr = split(/\|/,$nec_iprid);
                       }
                       else {
                               push(@iprscan_ipr,$nec_iprid);
                       }
               }
		foreach (@iprscan_ipr) {
			unless (exists ($iprscan_iprH{$_})) {
				push (@iprscan_iprA,$_);
				$iprscan_iprH{$_} += 1;
			}
		}

		# identify iprscan IPR_desc repeats and remove
		if ($ngodesc_iprdesc =~ /^#DESC#/) {
                       if ($niprdesc =~ m/\|/) {
                               @iprscan_iprdesc = split(/\|/,$niprdesc);
                       }
                       else {
                               push(@iprscan_iprdesc,$niprdesc);
                       }
               }
                foreach (@iprscan_iprdesc) {
                        unless (exists ($iprscan_iprdescH{$_})) {
                                push (@iprscan_iprdescA,$_);
                                $iprscan_iprdescH{$_} += 1;
                        }
                }

		# identify iprscan GO repeats and remove
		if ($ngotag =~ /^#GO#/) {
			if ($niprgo =~ m/\|/) {
				@iprscan_go= split(/\|/,$niprgo);
			}
			else {
				push (@iprscan_go,$niprgo);
			}	
		}	
		foreach (@iprscan_go) {
			unless(exists ($iprscan_goH{$_})) {
				push (@iprscan_goA,$_);
				$iprscan_goH{$_} += 1;
			}
		}
		# identify iprscan pathways repeats and remove
                if ($npathtag =~ /^#PATHWAY#/) {
                        if ($npathways =~ m/\|/) {
                                @iprscan_path= split(/\|/,$npathways);
                        }
                        else {
                                push (@iprscan_path,$npathways);
                        }
                }
                foreach (@iprscan_path) {
                        unless(exists ($iprscan_pathH{$_})) {
                                push (@iprscan_pathA,$_);
                                $iprscan_pathH{$_} += 1;
                        }
                }
                foreach (@iprscan_pathA) {  # extrat the EC number based on KEGG information
					if(/^KEGG:\s+\d+\+(.*)/) { # get only the EC:number from the line
						unless (exists ($iprscan_ecH{$1})) { # get only one copy of EC number
							push (@iprscan_ecA,"EC:$1");
							$iprscan_ecH{$1} += 1;
						}
					}
				}

		#print "$id\t$go_ipr\t$ec_iprid\t$godesc_iprdesc\t$nec_iprid\t$niprdesc\t$niprgo\t$npathways\n";
		local $"="|";

    # Was working as well, but without NULL for interproscan EC column field
		#print "$id\t@blast2go_goA\t$ec_iprid\t$godesc_iprdesc\t@iprscan_iprA\t@iprscan_iprdescA\t@iprscan_goA\t@iprscan_ecA\t$npathways\n";

    # Check if interproscan EC number is NULL or not
    print "$id\t@blast2go_goA\t$ec_iprid\t$godesc_iprdesc\t@iprscan_iprA\t@iprscan_iprdescA\t@iprscan_goA\t", scalar(@iprscan_ecA) == "" ? "NULL" : "@iprscan_ecA", "\t$npathways\n";


		local $"=" ";
		(@blast2go_go,@iprscan_ipr, @iprscan_go,@iprscan_iprdesc,@iprscan_path)=();
		(%blast2go_goH,%iprscan_iprH,%iprscan_goH,%iprscan_iprdescH,%iprscan_pathH,%iprscan_ecH)=();
		(@blast2go_goA,@iprscan_iprA,@iprscan_goA,@iprscan_iprdescA,@iprscan_pathA,@iprscan_ecA)=();
	}
	elsif ($id ne $nid) {
        	$boolean = 1;
        	#print "$id ne $nid\n";
		
		#print"$id\t@iprscan_iprA\t@blast2go_goA\t*mis*\n";	
		
		# identify blast2go GO repeats and remove
                if ($go_ipr =~ /^GO/) {
                       if ($go_ipr =~ m/\|/) {
                               @blast2go_go = split(/\|/,$go_ipr);
                       }
                       else {
                               push(@blast2go_go,$go_ipr);
                       }
                	foreach (@blast2go_go) {
                        	unless (exists ($blast2go_goH{$_})) {
                                	push (@blast2go_goA,$_);
                                	$blast2go_goH{$_} += 1;
                        	}
					}
					local $"="|";
                	print "$id\t@blast2go_goA\t$ec_iprid\t$godesc_iprdesc\tNULL\tNULL\tNULL\tNULL\tNULL\n";
					(@blast2go_go,@iprscan_ipr, @iprscan_go,@iprscan_iprdesc,@iprscan_path)=();
                    (%blast2go_goH,%iprscan_iprH,%iprscan_goH,%iprscan_iprdescH,%iprscan_pathH,%iprscan_ecH)=();
                    (@blast2go_goA,@iprscan_iprA,@iprscan_goA,@iprscan_iprdescA,@iprscan_pathA,@iprscan_ecA)=();

                }

                if ($go_ipr =~ /^#IPRID#/) {
			# identify iprscan IPR id repeats and remove
                       if ($ec_iprid =~ m/\|/) {
                               @iprscan_ipr = split(/\|/,$ec_iprid);
                       }
                       else {
                               push(@iprscan_ipr,$ec_iprid);
                       }

			# identify iprscan IPR_desc repeats and remove
                       if ($iprdesc =~ m/\|/) {
                               @iprscan_iprdesc = split(/\|/,$iprdesc);
                       }
                       else {
                               push(@iprscan_iprdesc,$iprdesc);
                       }

                	# identify iprscan GO repeats and remove
                        if ($iprgo =~ m/\|/) {
                                @iprscan_go= split(/\|/,$iprgo);
                        }
                        else {
                                push (@iprscan_go,$iprgo);
                        }
                
			# identify iprscan pathways repeats and remove
                        if ($pathways =~ m/\|/) {
                                @iprscan_path= split(/\|/,$pathways);
                        }
                        else {
                                push (@iprscan_path,$pathways);
                        }
                        foreach (@iprscan_ipr) {
                	        unless (exists ($iprscan_iprH{$_})) {
                        	        push (@iprscan_iprA,$_);
                                	$iprscan_iprH{$_} += 1;
                        	}
               		 }
                	foreach (@iprscan_iprdesc) {
                        	unless (exists ($iprscan_iprdescH{$_})) {
                                	push (@iprscan_iprdescA,$_);
                                	$iprscan_iprdescH{$_} += 1;
                        	}
                	}
                	foreach (@iprscan_go) {
                        	unless(exists ($iprscan_goH{$_})) {
                                	push (@iprscan_goA,$_);
                             	   $iprscan_goH{$_} += 1;
                        	}
                	}
          	      	foreach (@iprscan_path) {
                	        unless(exists ($iprscan_pathH{$_})) {
                        	        push (@iprscan_pathA,$_);
                                	$iprscan_pathH{$_} += 1;
                       	 	}
                	}
                	foreach (@iprscan_pathA) {  # extrat the EC number based on KEGG information
						if(/^KEGG:\s+\d+\+(.*)/) { # get only the EC:number from the line
							unless (exists ($iprscan_ecH{$1})) { # get only one copy of EC number
								push (@iprscan_ecA,"EC:$1");
								$iprscan_ecH{$1} += 1;
							}
						}
					}

					#print "$id\t$go_ipr\t$ec_iprid\t$godesc_iprdesc\t$iprdesc\t$gotag\t$iprgo\t$pathtag\t$pathways\t*mis*\n";
					local $"="|";
	                #print "$id\t@blast2go_goA\t$ec_iprid\t$godesc_iprdesc\t@iprscan_iprA\t@iprscan_iprdescA\t@iprscan_goA\t$pathways\n";
	                
                  # Was working as well, but without NULL for interproscan EC column field
                  #print "$id\tNULL\tNULL\tNULL\t@iprscan_iprA\t@iprscan_iprdescA\t@iprscan_goA\t@iprscan_ecA\t$pathways\n";

                  # Check if interproscan EC number is NULL or not
                  print "$id\tNULL\tNULL\tNULL\t@iprscan_iprA\t@iprscan_iprdescA\t@iprscan_goA\t", scalar(@iprscan_ecA) == "" ? "NULL" : "@iprscan_ecA", "\t$pathways\n";
                  
					#local $"=" ";
					(@blast2go_go,@iprscan_ipr, @iprscan_go,@iprscan_iprdesc,@iprscan_path)=();
					(%blast2go_goH,%iprscan_iprH,%iprscan_goH,%iprscan_iprdescH,%iprscan_pathH,%iprscan_ecH)=();
					(@blast2go_goA,@iprscan_iprA,@iprscan_goA,@iprscan_iprdescA,@iprscan_pathA,@iprscan_ecA)=();
        	}
	}
        else {
        	print "Problem with lines:\n$_\n$newline\n";
        }
}
close($inf);

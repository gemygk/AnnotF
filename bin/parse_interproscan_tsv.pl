#!/usr/bin/env perl


#parse the iprscan output file

# 	Gemy George Kaithakottil (gemy.kaithakottil@tgac.ac.uk)
# 	Version 0.1
#	Part of interproscan and blast2go pipeline

use strict;
#use warnings;

my $usage = "\nUsage: perl $0 <iprscan_file.tsv>

Note: The file needed to be sorted by first column (sort -k1,1 <filename>, will do)

";
my $file = shift or die $usage;

my $inf;
my $true=1;
my $for_pass="";
my ($id,$crc,$len,$method,$dbentry,$dbdesc,$start,$end,$eval,$status,$date,$ipr,$iprdesc,$go,$kegg);
my ($spos,$tpos,$coilspos,$phobiuspos,$signalpeukpos,$signalpgramppos,$signalpgramnpos); #for Seg and TMHMM regions
my (@pfamdbentry,@pfamipr,@pfamiprdesc,@pfamgo,@pfamkegg,
	@pirsfdbentry,@pirsfipr,@pirsfiprdesc,@pirsfgo,@pirsfkegg,
	@printsdbentry,@printsipr,@printsiprdesc,@printsgo,@printskegg,
	@profilescanipr,@profilescaniprdesc,@profilescango,
	@gene3ddbentry,@gene3dipr,@gene3diprdesc,@gene3dgo,@gene3dkegg,
	@pantherdbentry,@pantheripr,@pantheriprdesc,@panthergo,@pantherkegg,
	@segreg,
	@tmhmmreg,
	@superfamilydbentry,@superfamilyipr,@superfamilyiprdesc,@superfamilygo,@superfamilykegg,
	@hamapdbentry,@hamapipr,@hamapiprdesc,@hamapgo,@hamapkegg,
	@prositepatternsdbentry,@prositepatternsipr,@prositepatternsiprdesc,@prositepatternsgo,@prositepatternskegg,
	@prositeprofilesdbentry,@prositeprofilesipr,@prositeprofilesiprdesc,@prositeprofilesgo,@prositeprofileskegg,
	@prodomdbentry,@prodomipr,@prodomiprdesc,@prodomgo,@prodomkegg,
	@smartdbentry,@smartipr,@smartiprdesc,@smartgo,@smartkegg,
	@tigrfamdbentry,@tigrfamipr,@tigrfamiprdesc,@tigrfamgo,@tigrfamkegg,
	@patternscanipr,@patternscaniprdesc,@patternscango,
	@coilsreg,
	@phobiusdbentry,@phobiusreg,
	@signalpeukreg,@signalpeukdbentry,
	@signalpgrampreg,@signalpgrampdbentry,
	@signalpgramnreg,@signalpgramndbentry) =(); 

my (%gene3dh1,%gene3dh2,%gene3dh3,%gene3dh4,%gene3dh5,
	%pantherh1,%pantherh2,%pantherh3,%pantherh4,%pantherh5,
	%pfamh1,%pfamh2,%pfamh3,%pfamh4,%pfamh5,
	%pirsfh1,%pirsfh2,%pirsfh3,%pirsfh4,%pirsfh5,
	%printsh1,%printsh2,%printsh3,%printsh4,%printsh5,
	%smarth1,%smarth2,%smarth3,%smarth4,%smarth5,
	%superfamilyh1,%superfamilyh2,%superfamilyh3,%superfamilyh4,%superfamilyh5,
	%tigrfamh1,%tigrfamh2,%tigrfamh3,%tigrfamh4,%tigrfamh5,
	%hamaph1,%hamaph2,%hamaph3,%hamaph4,%hamaph5,
	%prositepatternsh1,%prositepatternsh2,%prositepatternsh3,%prositepatternsh4,%prositepatternsh5,
	%prositeprofilesh1,%prositeprofilesh2,%prositeprofilesh3,%prositeprofilesh4,%prositeprofilesh5,
	%prodomh1,%prodomh2,%prodomh3,%prodomh4,%prodomh5) =();
	
my (@mainid,@maindesc,@maingo,@mainkegg);
my (%mainidh,%maindesch,%maingoh,%mainkeggh);
my (@mainida,@maindesca,@maingoa,@mainkegga);

my @ids;
	
open $inf, "<" . $file;

while (<$inf>) {
	chomp();
	#ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/README.html
	##1####2###3#####4########5########6######7######8####9#####10######11###12#####13####14####15
	($id,$crc,$len,$method,$dbentry,$dbdesc,$start,$end,$eval,$status,$date,$ipr,$iprdesc,$go,$kegg) = split(/\t/);

	# create an if loop to skip newline
	
	# copy the first id
	if($true) {
		$for_pass=$id;
		$true=0;
	}
	if ($id ne $for_pass) {
		local $"="|";
		# here we try to print in the order iprids, desc, goterms, kegg
		#print"$for_pass\t#PANTHER#\t@pantherdbentry\t@pantheripr\t@pantheriprdesc\t@panthergo\t@pantherkegg\t#Pfam#\t@pfamdbentry\t@pfamipr\t@pfamiprdesc\t@pfamgo\t@pfamkegg\t#SUPERFAMILY#\t@superfamilydbentry\t@superfamilyipr\t@superfamilyiprdesc\t@superfamilygo\t@superfamilykegg\t#Gene3D#\t@gene3ddbentry\t@gene3dipr\t@gene3diprdesc\t@gene3dgo\t@gene3dkegg\t#PRINTS#\t@printsdbentry\t@printsipr\t@printsiprdesc\t@printsgo\t@printskegg\t#SMART#\t@smartdbentry\t@smartipr\t@smartiprdesc\t@smartgo\t@smartkegg\t#TIGRFAM#\t@tigrfamdbentry\t@tigrfamipr\t@tigrfamiprdesc\t@tigrfamgo\t@tigrfamkegg\t#PIRSF#\t@pirsfdbentry\t@pirsfipr\t@pirsfiprdesc\t@pirsfgo\t@pirsfkegg\t#ProSiteProfiles#\t@prositeprofilesdbentry\t@prositeprofilesipr\t@prositeprofilesiprdesc\t@prositeprofilesgo\t@prositeprofileskegg\t#ProSitePatterns#\t@prositepatternsdbentry\t@prositepatternsipr\t@prositepatternsiprdesc\t@prositepatternsgo\t@prositepatternskegg\t#Hamap#\t@hamapdbentry\t@hamapipr\t@hamapiprdesc\t@hamapgo\t@hamapkegg\t#ProDom#\t@prodomdbentry\t@prodomipr\t@prodomiprdesc\t@prodomgo\t@prodomkegg\t#Phobius#\t@phobiusdbentry\t@phobiusreg\t#TMHMM#\t@tmhmmreg\t#Coils#\t@coilsreg\t#SignalP_EUK#\t@signalpeukdbentry\t@signalpeukreg\t#SignalP_GRAM_POSITIVE#\t@signalpgrampdbentry\t@signalpgrampreg\t#SignalP_GRAM_NEGATIVE#\t@signalpgramndbentry\t@signalpgramnreg\t\n"; #main one
		
		
		@mainid = (@pantheripr,@pfamipr,@superfamilyipr,@gene3dipr,@printsipr,@smartipr,@tigrfamipr,@pirsfipr,@prositeprofilesipr,@prositepatternsipr,@hamapipr,@prodomipr);
		@maindesc = (@pantheriprdesc,@pfamiprdesc,@printsiprdesc,@superfamilyiprdesc,@gene3diprdesc,@smartiprdesc,@tigrfamiprdesc,@pirsfiprdesc,@prositeprofilesiprdesc,@prositepatternsiprdesc,@hamapiprdesc,@prodomiprdesc);
		@maingo = (@panthergo,@pfamgo,@superfamilygo,@gene3dgo,@printsgo,@smartgo,@tigrfamgo,@pirsfgo,@prositeprofilesgo,@prositepatternsgo,@hamapgo,@prodomgo);
		@mainkegg = (@pantherkegg,@pfamkegg,@superfamilykegg,@gene3dkegg,@printskegg,@smartkegg,@tigrfamkegg,@pirsfkegg,@prositeprofileskegg,@prositepatternskegg,@hamapkegg,@prodomkegg);
		
		foreach (@mainid) {
			unless (exists ($mainidh{$_})) {
				push (@mainida,$_);
				$mainidh{$_} += 1;
			}
		}
		foreach (@maindesc) {
			unless (exists ($maindesch{$_})) {
				push (@maindesca,$_);
				$maindesch{$_} += 1;
			}
		}
		foreach (@maingo) {
			unless (exists ($maingoh{$_})) {
				push (@maingoa,$_);
				$maingoh{$_} += 1;
			}
		}
		foreach (@mainkegg) {
			unless (exists ($mainkeggh{$_})) {
				push (@mainkegga,$_);
				$mainkeggh{$_} += 1;
			}
		}
		
		#copy from here to eof
		if($for_pass =~ /\|/) { 																				# checking if contigs has piping
			@ids = split(/\|/,$for_pass); 																		# store the split contigs to array
			foreach (@ids) {																					# call each contigs straight away
				#print"$_\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n";		# print the contigs and the rest
				print"$_\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n" unless ((scalar(@mainida) == "") && (scalar(@maindesca) == "") && (scalar(@maingoa) == "") && (scalar(@mainkegga) == "")); #main one
			}
			@ids=();																							# empty the array
		}
		else {																									# normal uniq contigs printed straight away
			#print"$for_pass\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n"; #main one
			print"$for_pass\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n" unless ((scalar(@mainida) == "") && (scalar(@maindesca) == "") && (scalar(@maingoa) == "") && (scalar(@mainkegga) == "")); #main one
		}
		
		#empty variables created for the final output
		(@mainid,@maindesc,@maingo,@mainkegg) =();
		(%mainidh,%maindesch,%maingoh,%mainkeggh)=();
		(@mainida,@maindesca,@maingoa,@mainkegga)=();
		
		$for_pass=$id;
		#empty arrays
		(@printsdbentry,@printsipr,@printsiprdesc,@printsgo,@printskegg,
		@pfamdbentry,@pfamipr,@pfamiprdesc,@pfamgo,@pfamkegg,
		@pirsfdbentry,@pirsfipr,@pirsfiprdesc,@pirsfgo,@pirsfkegg,
		@profilescanipr,@profilescaniprdesc,@profilescango,
		@gene3ddbentry,@gene3dipr,@gene3diprdesc,@gene3dgo,@gene3dkegg,
		@pantherdbentry,@pantheripr,@pantheriprdesc,@panthergo,@pantherkegg,
		@segreg,
		@tmhmmreg,
		@superfamilydbentry,@superfamilyipr,@superfamilyiprdesc,@superfamilygo,@superfamilykegg,
		@hamapdbentry,@hamapipr,@hamapiprdesc,@hamapgo,@hamapkegg,
		@prositepatternsdbentry,@prositepatternsipr,@prositepatternsiprdesc,@prositepatternsgo,@prositepatternskegg,
		@prositeprofilesdbentry,@prositeprofilesipr,@prositeprofilesiprdesc,@prositeprofilesgo,@prositeprofileskegg,
		@prodomdbentry,@prodomipr,@prodomiprdesc,@prodomgo,@prodomkegg,
		@smartdbentry,@smartipr,@smartiprdesc,@smartgo,@smartkegg,
		@tigrfamdbentry,@tigrfamipr,@tigrfamiprdesc,@tigrfamgo,@tigrfamkegg,
		@patternscanipr,@patternscaniprdesc,@patternscango,
		@coilsreg,
		@phobiusdbentry,@phobiusreg,
		@signalpeukreg,@signalpeukdbentry,
		@signalpgrampreg,@signalpgrampdbentry,
		@signalpgramnreg,@signalpgramndbentry) =(); 
		#empty hashes
		(%gene3dh1,%gene3dh2,%gene3dh3,%gene3dh4,%gene3dh5,
		%pantherh1,%pantherh2,%pantherh3,%pantherh4,%pantherh5,
		%pfamh1,%pfamh2,%pfamh3,%pfamh4,%pfamh5,
		%pirsfh1,%pirsfh2,%pirsfh3,%pirsfh4,%pirsfh5,
		%printsh1,%printsh2,%printsh3,%printsh4,%printsh5,
		%smarth1,%smarth2,%smarth3,%smarth4,%smarth5,
		%superfamilyh1,%superfamilyh2,%superfamilyh3,%superfamilyh4,%superfamilyh5,
		%tigrfamh1,%tigrfamh2,%tigrfamh3,%tigrfamh4,%tigrfamh5,
		%hamaph1,%hamaph2,%hamaph3,%hamaph4,%hamaph5,
		%prositepatternsh1,%prositepatternsh2,%prositepatternsh3,%prositepatternsh4,%prositepatternsh5,
		%prositeprofilesh1,%prositeprofilesh2,%prositeprofilesh3,%prositeprofilesh4,%prositeprofilesh5,
		%prodomh1,%prodomh2,%prodomh3,%prodomh4,%prodomh5) =();
		#empty scalars
		($spos,$tpos,$coilspos,$phobiuspos,$signalpeukpos,$signalpgramppos,$signalpgramnpos)="";
	}
	if($id eq $for_pass) {
		#print "$id,$crc,$len,$method,$dbentry,$dbdesc,$start,$end,$eval,$status,$date,$ipr,$iprdesc,$go\n"; #testing
		
		if ($method eq "Gene3D") { #1#
			if($dbentry ne "") {
				unless(exists($gene3dh1{$dbentry})) {
					push(@gene3ddbentry,$dbentry); # get database members entry
					 $gene3dh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($gene3dh2{$ipr})) {
					push(@gene3dipr,$ipr);
					$gene3dh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($gene3dh3{$iprdesc})) {
					push(@gene3diprdesc,$iprdesc);
					$gene3dh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($gene3dh4{$go})) {
					push(@gene3dgo,$go);
					$gene3dh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($gene3dh5{$kegg})) {
					push(@gene3dkegg,$kegg);
					$gene3dh5{$kegg} +=1;
				}
			}
			#@gene3dgo,@gene3dkegg
			
		}
		if ($method eq "Coils") { #2#
			$coilspos="$start-$end";
			push(@coilsreg,$coilspos); # get start and stop of domain match
		}
		if ($method eq "PANTHER") { #3#
			if($dbentry ne "") {
				unless(exists($pantherh1{$dbentry})) {
					push(@pantherdbentry,$dbentry); # get database members entry
					$pantherh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($pantherh2{$ipr})) {
					push(@pantheripr,$ipr);
					$pantherh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($pantherh3{$iprdesc})) {
					push(@pantheriprdesc,$iprdesc);
					$pantherh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($pantherh4{$go})) {
					push(@panthergo,$go);
					$pantherh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($pantherh5{$kegg})) {
					push(@pantherkegg,$kegg);
					$pantherh5{$kegg} +=1;
				}
			}
		}
		if ($method eq "Pfam") {  #4#
			if($dbentry ne "") {
				unless(exists($pfamh1{$dbentry})) {
					push(@pfamdbentry,$dbentry);
					$pfamh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($pfamh2{$ipr})) {
					push(@pfamipr,$ipr);
					$pfamh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($pfamh3{$iprdesc})) {
					push(@pfamiprdesc,$iprdesc);
					$pfamh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($pfamh4{$go})) {
					push(@pfamgo,$go);
					$pfamh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($pfamh5{$kegg})) {
					push(@pfamkegg,$kegg);
					$pfamh5{$kegg} +=1;
				}
			}
		}
		if ($method eq "Phobius") { #5#
			$phobiuspos="$start-$end";
			push(@phobiusreg,$phobiuspos); # get start and stop of domain match
			push(@phobiusdbentry,$dbentry); # get database members entry
			#@phobiusdbentry,@phobiusreg
		}
		if ($method eq "PIRSF") {  #6#
			if($dbentry ne "") {
				unless(exists($pirsfh1{$dbentry})) {
					push(@pirsfdbentry,$dbentry);
					$pirsfh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($pirsfh2{$ipr})) {
					push(@pirsfipr,$ipr);
					$pirsfh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($pirsfh3{$iprdesc})) {
					push(@pirsfiprdesc,$iprdesc);
					$pirsfh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($pirsfh4{$go})) {
					push(@pirsfgo,$go);
					$pirsfh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($pirsfh5{$kegg})) {
					push(@pirsfkegg,$kegg);
					$pirsfh5{$kegg} +=1;
				}
			}
		}
		if($method eq "PRINTS")	{    #7#
			#print "$id\t$ipr\t$iprdesc\t$go\n";
			#push(@printsipr,$ipr); # get InterPro entry
			#push(@printsiprdesc,$iprdesc); #get InterPro entry description
			#push(@printsgo,$go); # get GO (gene ontology) description
			if($dbentry ne "") {
				unless(exists($printsh1{$dbentry})) {
					push(@printsdbentry,$dbentry);
					$printsh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($printsh2{$ipr})) {
					push(@printsipr,$ipr);
					$printsh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($printsh3{$iprdesc})) {
					push(@printsiprdesc,$iprdesc);
					$printsh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($printsh4{$go})) {
					push(@printsgo,$go);
					$printsh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($printsh5{$kegg})) {
					push(@printskegg,$kegg);
					$printsh5{$kegg} +=1;
				}
			}
		}
		if ($method eq "SignalP_EUK") { #8#
			$signalpeukpos="$start-$end";
			push(@signalpeukreg,$signalpeukpos); # get start and stop of domain match
			push(@signalpeukdbentry,$dbentry); # get database members entry
			#@signalpeukreg,@signalpeukdbentry
		}
		if ($method eq "SignalP_GRAM_POSITIVE") { #8#
			$signalpgramppos="$start-$end";
			push(@signalpgrampreg,$signalpgramppos); # get start and stop of domain match
			push(@signalpgrampdbentry,$dbentry); # get database members entry
			#@signalpgrampreg,@signalpgrampdbentry
		}
		if ($method eq "SignalP_GRAM_NEGATIVE") { #8#
			$signalpgramnpos="$start-$end";
			push(@signalpgramnreg,$signalpgramnpos); # get start and stop of domain match
			push(@signalpgramndbentry,$dbentry); # get database members entry
			#@signalpgramnreg,@signalpgramndbentry
		}
		if($method eq "SMART")	{    #9#
			if($dbentry ne "") {
				unless(exists($smarth1{$dbentry})) {
					push(@smartdbentry,$dbentry);
					$smarth1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($smarth2{$ipr})) {
					push(@smartipr,$ipr);
					$smarth2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($smarth3{$iprdesc})) {
					push(@smartiprdesc,$iprdesc);
					$smarth3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($smarth4{$go})) {
					push(@smartgo,$go);
					$smarth4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($smarth5{$kegg})) {
					push(@smartkegg,$kegg);
					$smarth5{$kegg} +=1;
				}
			}
		}
		if($method eq "SUPERFAMILY")	{    #10#
			if($dbentry ne "") {
				unless(exists($superfamilyh1{$dbentry})) {
					push(@superfamilydbentry,$dbentry);
					$superfamilyh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($superfamilyh2{$ipr})) {
					push(@superfamilyipr,$ipr);
					$superfamilyh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($superfamilyh3{$iprdesc})) {
					push(@superfamilyiprdesc,$iprdesc);
					$superfamilyh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($superfamilyh4{$go})) {
					push(@superfamilygo,$go);
					$superfamilyh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($superfamilyh5{$kegg})) {
					push(@superfamilykegg,$kegg);
					$superfamilyh5{$kegg} +=1;
				}
			}
		}
		if($method eq "TIGRFAM")	{    #11#
			if($dbentry ne "") {
				unless(exists($tigrfamh1{$dbentry})) {
					push(@tigrfamdbentry,$dbentry);
					$tigrfamh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($tigrfamh2{$ipr})) {
					push(@tigrfamipr,$ipr);
					$tigrfamh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($tigrfamh3{$iprdesc})) {
					push(@tigrfamiprdesc,$iprdesc);
					$tigrfamh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($tigrfamh4{$go})) {
					push(@tigrfamgo,$go);
					$tigrfamh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($tigrfamh5{$kegg})) {
					push(@tigrfamkegg,$kegg);
					$tigrfamh5{$kegg} +=1;
				}
			}
		}
		if ($method eq "TMHMM") { #12#
			$tpos="$start-$end";
			push(@tmhmmreg,$tpos); # get start and stop of domain match
		}
		
		if($method eq "Hamap")	{    #13#
			if($dbentry ne "") {
				unless(exists($hamaph1{$dbentry})) {
					push(@hamapdbentry,$dbentry);
					$hamaph1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($hamaph2{$ipr})) {
					push(@hamapipr,$ipr);
					$hamaph2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($hamaph3{$iprdesc})) {
					push(@hamapiprdesc,$iprdesc);
					$hamaph3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($hamaph4{$go})) {
					push(@hamapgo,$go);
					$hamaph4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($hamaph5{$kegg})) {
					push(@hamapkegg,$kegg);
					$hamaph5{$kegg} +=1;
				}
			}
		}
		if($method eq "ProDom")	{    #14#
			if($dbentry ne "") {
				unless(exists($prodomh1{$dbentry})) {
					push(@prodomdbentry,$dbentry);
					$prodomh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($prodomh2{$ipr})) {
					push(@prodomipr,$ipr);
					$prodomh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($prodomh3{$iprdesc})) {
					push(@prodomiprdesc,$iprdesc);
					$prodomh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($prodomh4{$go})) {
					push(@prodomgo,$go);
					$prodomh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($prodomh5{$kegg})) {
					push(@prodomkegg,$kegg);
					$prodomh5{$kegg} +=1;
				}
			}
		}
		if($method eq "ProSitePatterns")	{    #15#
			if($dbentry ne "") {
				unless(exists($prositepatternsh1{$dbentry})) {
					push(@prositepatternsdbentry,$dbentry);
					$prositepatternsh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($prositepatternsh2{$ipr})) {
					push(@prositepatternsipr,$ipr);
					$prositepatternsh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($prositepatternsh3{$iprdesc})) {
					push(@prositepatternsiprdesc,$iprdesc);
					$prositepatternsh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($prositepatternsh4{$go})) {
					push(@prositepatternsgo,$go);
					$prositepatternsh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($prositepatternsh5{$kegg})) {
					push(@prositepatternskegg,$kegg);
					$prositepatternsh5{$kegg} +=1;
				}
			}
		}
		if($method eq "ProSiteProfiles")	{    #16#
			if($dbentry ne "") {
				unless(exists($prositeprofilesh1{$dbentry})) {
					push(@prositeprofilesdbentry,$dbentry);
					$prositeprofilesh1{$dbentry} +=1;
				}
			}
			if($ipr ne "") {
				unless(exists($prositeprofilesh2{$ipr})) {
					push(@prositeprofilesipr,$ipr);
					$prositeprofilesh2{$ipr} +=1;
				}
			}
			if($iprdesc ne "") {
				unless(exists($prositeprofilesh3{$iprdesc})) {
					push(@prositeprofilesiprdesc,$iprdesc);
					$prositeprofilesh3{$iprdesc} +=1;
				}
			}
			if($go ne "") {
				unless(exists($prositeprofilesh4{$go})) {
					push(@prositeprofilesgo,$go);
					$prositeprofilesh4{$go} +=1;
				}
			}
			if($kegg ne "") {
				unless(exists($prositeprofilesh5{$kegg})) {
					push(@prositeprofileskegg,$kegg);
					$prositeprofilesh5{$kegg} +=1;
				}
			}
		}
	}
	if(eof($inf)) {
		local $"="|";
		foreach (@mainid) {
			unless (exists ($mainidh{$_})) {
				push (@mainida,$_);
				$mainidh{$_} += 1;
			}
		}
		foreach (@maindesc) {
			unless (exists ($maindesch{$_})) {
				push (@maindesca,$_);
				$maindesch{$_} += 1;
			}
		}
		foreach (@maingo) {
			unless (exists ($maingoh{$_})) {
				push (@maingoa,$_);
				$maingoh{$_} += 1;
			}
		}
		foreach (@mainkegg) {
			unless (exists ($mainkeggh{$_})) {
				push (@mainkegga,$_);
				$mainkeggh{$_} += 1;
			}
		}
		
		#copy from here to eof
		if($for_pass =~ /\|/) { 																				# checking if contigs has piping
			@ids = split(/\|/,$for_pass); 																		# store the split contigs to array
			foreach (@ids) {																					# call each contigs straight away
				#print"$_\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n";		# print the contigs and the rest
				print"$_\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n" unless ((scalar(@mainida) == "") && (scalar(@maindesca) == "") && (scalar(@maingoa) == "") && (scalar(@mainkegga) == "")); #main one
			}
			@ids=();																							# empty the array
		}
		else {																									# normal uniq contigs printed straight away
			#print"$for_pass\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n"; #main one
			print"$for_pass\t#IPRID#\t@mainida\t#DESC#\t@maindesca\t#GO#\t@maingoa\t#PATHWAY#\t@mainkegga\n" unless ((scalar(@mainida) == "") && (scalar(@maindesca) == "") && (scalar(@maingoa) == "") && (scalar(@mainkegga) == "")); #main one
		}
	}
}
close($inf);

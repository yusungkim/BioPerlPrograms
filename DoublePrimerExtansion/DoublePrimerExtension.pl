####### DoublePrimerExtension.pl
####### 15.Sep.2011 YuSung Kim
####### usage: DoublePrimerExtension.pl [options] <Input_File> [-o Output_File] [-e Enzyme]
####### 하는일 : 입력된 풀 시퀀스 DNA(N:ACGT, S:CG 등의 문자도 포함, IUB ambiguity codes)를 합성할 수있는 두개의 프라이머를 생성한다.
####### 알고리즘 : 프라이머 Fw와 Rv의 3말단의 일부를 상보적으로 만든 후, 나머지부분은 랜덤으로 생성한다.
#######          그 다음 랜덤으로 생성된 프라이머들을 alignment하여 그 값들을 비교한다.
#######          3말단의 상보적 배열이 가장 높은 스코어로 결합되는 페어들을 찾아 낸후,
#######          그 중에서 일어날 수 있는 다른 부분에서의 결합자유에너지가 가장높은(결합잘 일어나지않는) 페어를 찾아낸다.
#######          이것을 통해 3말단의 상보적배열의 결합이 특이적으로 일어나는 프라이머 페어를 생성한다.

sub version {
  print "\n";
  print "###########################################################\n";
  print "##         DoublePrimerExtension Version 1.1             ##\n";
  print "##     .............................................     ##\n";
  print "##             All Copywrights are reserved by           ##\n";
  print "##           Yusung KIM, 2011.09.14 ~ 2011.09.23         ##\n";
  print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com  ##\n";
  print "###########################################################\n";
}
sub usage{
    print "\nusage: perl $0 [options] <Input File> [Output File]\n";
    print "ex)\n perl DoublePrimerExtension.pl -v DNAseq.txt -o result.txt -e KpnI -e KpnI\n";
    print "options:\n";
    print "	-v	-> display the process\n";
    print "	-d	-> display info for debug\n";
    print "	-V	-> display version\n";
    print "	-o	-> specify output file to save the results\n";
    print "		   ex) -o output.txt\n";
    print "	-e	-> specify restruction enzyme not to cut by\n";
    print "		   ex) -e KpnI\n";
    print "	-h/?	-> display this\n";
}


################################### MAIN #################################

#!/usr/bin/perl -w
use warnings;
use strict;

# set flags as default
my (%flags) = (
         "visual" => "false", ## donot display the process.
         "debug"  => "false", ## donot display info for debug.
         "input"  => "false", # input_file existence
         "output" => "false",	# output_file existence
         "restruction"	=> "true"	# specify restruction enzyme not to cut
	       );
	       
my (%enzymes); #DB read from file, RestructionEnzymeList, ex) "KpnI" => "GGTACC"
my @enzymesNotToCut;

######## DATA
# Default Hydrophilic Sequence
my ($fullPhilicSeq) = "GT TGG CAC GTC ATA TGC TCT TCA GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTN CCN GGN CGN GGN GTN CCN GGN GTN GGN GTA GGA AGA GTA ACG AGC TCA TGC";
# Default Hydrophobic Sequence
my ($fullPhobicSeq) = "GTT GGC ACG TGA GCT CCT CTT CA ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN GTN CCN GGN TAY GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN GTN CCN GGN TAY GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN GTN CCN GGN TAY GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN GTN CCN GGN TAY GGN ATH CCN GGN GTN GGN ATH CCN GGN GTN GGN ATH GGA AGA GTA ACG GTA CCA TGC";
my $input_file = '';
my $output_file = 'result.txt'; #default output file

######## parameters
#my $isPhilic = 1;
my ($sizeOfOverLap) = 25;	# OverLap Size, default 25 bp
my ($numOfFw) = 5;			# Number of Random Forward primers
my ($numOfRv) = 5;			# Reverse

######## temporalary data storages
my ($fullSequence);				# for protype of full seq.
my ($templateFw, $templateRv);	# for protype of primer seq.
my ($randFwSeq, $randRvSeq);	# for random generation of seq.
my (@ForwardPrimers);			# list of primers (1 demention)
my (@ReversePrimers);			# list of primers (2 demention) !!
my (@Candidates);				# for selecting primer pairs
my (@partialFwScore, @partialRvScore);	# for scoring at 2nd selection
my (@firstScoreBoard);

my ($sizeOfFullSeq) = 0;
my ($sizeOfPrimerFw, $sizeOfPrimerRv);
my ($ii, $jj);
#my ($randomNucleotide) = rand_nucleotide_N();
my $leastScore = 0;
my $leastScoredFw = '';
my $leastScoredRv = '';

####### CHECKING OPTIONS

if( @ARGV == 0 ) {
	version();
	usage();
	exit;
}

if( defined($ARGV[0]) && ($ARGV[0] =~ /^-\S*V\S*/) ) {
    version();
    exit;
}
if ( defined($ARGV[0]) && ($ARGV[0] =~ /^-\S*[h?]\S*/) ) {
	usage();
	exit;
}

my $temp = '';
while( defined($temp = shift) ) {
    
    # -v : visual option
	if($temp =~ /^-\S*v\S*/){
		$flags{"visual"} = "true";
	}
	# -d : debug option
	if($temp =~ /^-\S*d\S*/){
		$flags{"debug"} = "true";
	}
	# -o : output file
	if($temp =~ /^-\S*o\S*/){
		$flags{"output"} = "true";
		$output_file = shift;
		open(OUTPUT_FILE, ">", $output_file) || die "Cannot open:\n\t$!";
	}
	# -e : restruction enzyme
	if($temp =~ /^-\S*e\S*/){
		$flags{"restruction"} = "true";
		my $enzymeName = shift;
		# check if the enzname fllows.
		if ((!defined $enzymeName) || ($enzymeName =~ /^-/)) {
			print "Wrong option usage:\n";
			print "\tEnzyme name should be followed after -e option\n";
			exit;
		}
		push @enzymesNotToCut, $enzymeName;
	}
	
	# Wrong options
	if($temp =~ /^-[vdoe]*([^dvoe])+[dvoe]*$/){
		print "Wrong option used: $1\n";
		usage();
		exit;
	}
	
	####### CHECKING INPUT FILES
	# read seq. from source file
	if( defined($temp) && $temp =~ /^[^-]/){
		
		$input_file = $temp;
    	$flags{'input'} = 'true';
    	open(INPUT_FILE, $input_file) or die "Cannot open:\n\t$!";
    	
    	my $contents = '';
    	while (my $line = <INPUT_FILE>) {
    		$contents .= $line;
    	}
    	$fullSequence = $contents;
    	
    	close(INPUT_FILE);		
    }
}


####################### read restruction enzyme from file
open(RES, "RestructionEnzymeList.txt") || die "Cannot open:\n\t$!";
	while ( <RES> ) {
		chomp;
		$_ =~ /^(\S+)\s*,\s*(\S+)\s*/;
		my $name = $1;
		my $seq = $2;
		$seq =~ s/[^a-zA-Z]//;
		$enzymes{$name} = $seq;
	}
close(RES);
my $num = %enzymes;
print "$num Restruction Enzymes were read.\n";
		
####################### pre-process
# remove everything excetp DNA
$fullSequence =~ s/[^AGCTNYRSWHVDB]//ig;

# remove spaces
#$fullSequence =~ s/[\s]//g;

# remove numbers
#$fullSequence =~ s/\d//g;

# remove non-alphabet
#$fullSequence =~ s/[^A-Za-z]//g;

# change to Upper case
$fullSequence = uc $fullSequence; # $fullSequence =~ tr/[a-z]/[A-Z]/;


####################### set Primers
$sizeOfFullSeq = length $fullSequence;
if( ($sizeOfFullSeq + $sizeOfOverLap) %2 == 0) {
	$sizeOfPrimerFw = ($sizeOfFullSeq + $sizeOfOverLap) / 2;
	$sizeOfPrimerRv = $sizeOfPrimerFw;
}else{
	$sizeOfPrimerFw = int (($sizeOfFullSeq + $sizeOfOverLap) / 2);
	$sizeOfPrimerRv = $sizeOfPrimerFw + 1;
}

$templateFw = substr $fullSequence, 0, $sizeOfPrimerFw;
$templateRv = substr $fullSequence, -$sizeOfPrimerRv;


####################### Random primer generation
srand(time);
for (my $ctr1 = 0; $ctr1 < $numOfFw; $ctr1++) {
	
	# synthesis front part
	my $dummy1 = substr ($templateFw, 0, -$sizeOfOverLap);
	$dummy1 =~ s/N/&rand_nucleotide_weakN/ge; # weak seq for the front part
	$dummy1 =~ s/R/&rand_nucleotide_weakR/ge; # weak seq for the front part, R : A, G
	$dummy1 =~ s/H/&rand_nucleotide_weakH/ge; # weak seq for H: T, C, A
	$dummy1 =~ s/Y/&rand_nucleotide_weakY/ge; # Y : C, T
	$dummy1 =~ s/S/&rand_nucleotide_S/ge;
	$dummy1 =~ s/W/&rand_nucleotide_W/ge;
	$dummy1 =~ s/V/&rand_nucleotide_weakV/ge;
	$dummy1 =~ s/D/&rand_nucleotide_weakD/ge;
	$dummy1 =~ s/B/&rand_nucleotide_weakB/ge;
	
	# synthesis overlap part
	my $dummy2 = substr ($templateFw, -$sizeOfOverLap);
	$dummy2 =~ s/N/&rand_nucleotide_strongN/ge; # strong seq for the overlap part
	$dummy2 =~ s/R/&rand_nucleotide_strongR/ge; # weak seq for the front part, R : A, G
	$dummy2 =~ s/H/&rand_nucleotide_strongH/ge; # weak seq for H: T, C, A
	$dummy2 =~ s/Y/&rand_nucleotide_strongY/ge; # Y : C, T
	$dummy2 =~ s/S/&rand_nucleotide_S/ge;
	$dummy2 =~ s/W/&rand_nucleotide_W/ge;
	$dummy2 =~ s/V/&rand_nucleotide_strongV/ge;
	$dummy2 =~ s/D/&rand_nucleotide_strongD/ge;
	$dummy2 =~ s/B/&rand_nucleotide_strongB/ge;
	
	# synthesis full forward seq
	$randFwSeq = $dummy1.$dummy2;
	
	
	# check for Restruction Enzyme
	my $isCutSiteFound = "false";
	
	if ($flags{'restruction'} eq 'true') {
		# Restruction Enzyme site deletion
		foreach my $name (@enzymesNotToCut) {
			my $pattern = $enzymes{$name};
			if( defined $pattern ) {
				if( $randFwSeq =~ /$pattern/ ) {
					$ctr1--;
					$isCutSiteFound = "true";
					print "======================= $pattern for $name are found! -> Deleted!\n";
				}
			}
		}
		
		# if there was restruction enzyme site then do it again.
		if ($isCutSiteFound eq "true") {
			next;
		}
	}
	
	# save in primer list
	$ForwardPrimers[$ctr1] = $randFwSeq;
	
	# find reverse primer
	for (my $ctr2 = 0; $ctr2 < $numOfRv; $ctr2++) {
		# copy the overlap part from forward primer
		$randRvSeq = substr($ForwardPrimers[$ctr1], -$sizeOfOverLap);
		
		# set the last part of protype of reverse primer
		$randRvSeq .= substr($templateRv, $sizeOfOverLap);
		
		# A or T for remained seq.
		$randRvSeq =~ s/N/&rand_nucleotide_weakN/ge; # weak seq for the front part
		$randRvSeq =~ s/R/&rand_nucleotide_weakR/ge; # weak seq for the front part, R : A, G
		$randRvSeq =~ s/H/&rand_nucleotide_weakH/ge; # weak seq for H: T, C, A
		$randRvSeq =~ s/Y/&rand_nucleotide_weakY/ge; # Y : C, T
		$randRvSeq =~ s/S/&rand_nucleotide_S/ge;
		$randRvSeq =~ s/W/&rand_nucleotide_W/ge;
		$randRvSeq =~ s/V/&rand_nucleotide_weakV/ge;
		$randRvSeq =~ s/D/&rand_nucleotide_weakD/ge;
		$randRvSeq =~ s/B/&rand_nucleotide_weakB/ge;
		
		# check for Restruction Enzyme
		if ($flags{'restruction'} eq 'true') {
			my $isCutSiteFound = "false";
			
			# Restruction Enzyme site deletion
			foreach my $name (@enzymesNotToCut) {
				my $pattern = $enzymes{$name};
				if( defined $pattern ) {
					if( $randRvSeq =~ /$pattern/ ) {
						$ctr2--;
						$isCutSiteFound = "true";
						print "======================= $pattern for $name are found! -> Deleted!\n";
					}
				}
			}
			
			if ($isCutSiteFound eq "true") {
				next;
			}
		}
		
		# save in primer list
		$ReversePrimers[$ctr1][$ctr2] = $randRvSeq;
		
		# to complementary seq
		#$randRvSeq =~ tr/ATGCatgc/TACGtacg/;
		
		# to reverse seq
		#$randRvSeq = reverse $randRvSeq;
	
		# Alignment
		if ( $flags{"visual"} eq "true" ) {
			PrintOrDisp ("---------------------------------------------------------------------\n");
		}
		
		# Find Candidate (First Screening)
		# Only OverLap part aligmnented
		my $score = getScore ($ForwardPrimers[$ctr1], $ReversePrimers[$ctr1][$ctr2], $ctr1, $ctr2);
		$firstScoreBoard[$ctr1][$ctr2] = $score;
		my $ii = $ctr1 + 1;
		my $jj = $ctr2 + 1;
		
		PrintOrDisp ("First Selection => Fw$ii x Rv$jj ");
		if( $Candidates[$ctr1][$ctr2] == 1 ) {	## @Candidates was determined in sub function.
			PrintOrDisp ( "-> OK, Score:$score\n" );
		}else{
			PrintOrDisp ( "-> Not Good, Score:$score\n" );
		}
		if( $flags{'visual'} eq 'true' ) {
			my $buf = $ReversePrimers[$ctr1][$ctr2];
			$buf = reverse $buf;
			$buf =~ tr/ATGCatgc/TACGtacg/;
			PrintOrDisp ( "Fw$ii:\t$ForwardPrimers[$ctr1]\nRv$ii-$jj:\t$buf\n\n" );
		}
		
		
		# Find Candidate (Second Screening) within 1st screened pair
		if( $Candidates[$ctr1][$ctr2] == 1 ) {
			my $offset = int($sizeOfOverLap/2);
			
			# partial Fw vs Rv
			my $partialSeq = substr ($ForwardPrimers[$ctr1], 0, -$offset);
			my $score = getScore ($partialSeq, $ReversePrimers[$ctr1][$ctr2]);
			$partialFwScore[$ctr1][$ctr2] = $score;
			PrintOrDisp ( "\n" );
			PrintOrDisp ( "   Second Selection => partialFw x Rv, Score:$score\n");
			
			# Fw vs partial Rv
			$partialSeq = substr ($ReversePrimers[$ctr1][$ctr2], $offset);
			$score = getScore ($randFwSeq, $partialSeq);
			$partialRvScore[$ctr1][$ctr2] = $score;
			PrintOrDisp ("   Second Selection => Fw x partialRv, Socre:$score\n\n");
		}
	}
}

# Statics
my @GoodCandidates = ();
my @SortByScore1 = ();
my @SortByScore2 = ();
my $score = 0;
my $ID = 0;

# Gather into new array
#print " Finding most appropriate seq pair.......";
$ID = 0;
for(my $i=0; $i<$numOfFw; $i++){
	for (my $j=0; $j<$numOfRv; $j++) {
 		if($Candidates[$i][$j] == 1) {
 			push (@GoodCandidates, "$i\t$j"); #Offset of Fw & Rv
 			
 			$score = $partialFwScore[$i][$j];
 			my $formattedScore = sprintf ("%08d", $score);
 			push (@SortByScore1, "$formattedScore\t$ID\t$i\t$j"); #ID:as order of gathered pair, Fw, Rv
 			
 			$score = $partialRvScore[$i][$j];
 			$formattedScore = sprintf ("%08d", $score);
 			push (@SortByScore2, "$formattedScore\t$ID\t$i\t$j");
 			
 			$ID++;
 		}
 	}
}

if ($ID == 0) {
	PrintOrDisp( "There is no good candidate\n" );
	PrintOrDisp( "Try once more\n");
	PrintOrDisp( "My time is both limited and valuable\n");
	PrintOrDisp("In romance, as in show business, always leave them wanting more\n");
	exit;
}

# Sort by each soring system
@SortByScore1 = sort {$a cmp $b} @SortByScore1; #numeric sort
@SortByScore2 = sort {$a cmp $b} @SortByScore2; #numeric sort

# Sort by orders
my @totalOrder = ();
for (my $i=0; $i < @GoodCandidates; $i++) {
	$ID = $i;
	
	#find id
	my $info;
	my $order;
	for (my $j=0; $j < @SortByScore1; $j++) {
		$info = $SortByScore1[$j];
		$info =~ /(\d+)\t(\d+)\t(\d+)\t(\d+)/; # Score, ID, Fw, Rv
		if($2 == $ID) {
			$order = $j;
			$score = $1;
		}
	}
	for (my $j=0; $j < @SortByScore2; $j++) {
		$info = $SortByScore2[$j];
		$info =~ /(\d+)\t(\d+)\t(\d+)\t(\d+)/; # Score, ID, Fw, Rv
		if($2 == $ID) {
			$order += $j;
			$score += $1;
		}
	}
	my $formattedString = sprintf("%05d", $order);
	push (@totalOrder, "$formattedString\t$ID\t$score"); #order, ID, score
}
@totalOrder = sort {$a cmp $b} @totalOrder; #alphabetical sort


#################################### Final Report
# display all sequences if visual option is on
if($flags{'visual'} eq 'true') {
	PrintOrDisp ("\nSequences...\n");
	for(my $i=0; $i<@ForwardPrimers; $i++) {
		my $I = $i+1;
		PrintOrDisp ("Fw$I $ForwardPrimers[$i]\n");
	
		for(my $j=0; $j<$numOfRv; $j++) {
			my $J = $j+1;
			my $primerRvSeq = $ReversePrimers[$i][$j];
	
			# reverse complementary seq
			$primerRvSeq = reverse $primerRvSeq;
			$primerRvSeq =~ tr/ATGCatgc/TACGtacg/;
	
			PrintOrDisp (" Rv$I-$J $primerRvSeq\n");
		}
	}
}

my $comb = $numOfFw * $numOfRv;
my $numOfGoodCandidates = @GoodCandidates;

# Get Best the ID of primer pair
$totalOrder[0] =~ /(\d+)\t(\d+)\t(\d+)/; # order, ID, score
$ID = $2;

# Get Best primer seq
$GoodCandidates[$ID] =~ /(\d+)\t(\d+)/;	# offset of forward/reverse primer
my $primerFwSeq = $ForwardPrimers[$1];
my $primerRvSeq = $ReversePrimers[$1][$2];
$ii = $1+1;
$jj = $2+1;

PrintOrDisp( "=====================================================================\n" );
PrintOrDisp( "<Summary>\n" );
PrintOrDisp( "We compaired $numOfFw random Forward sequences\n");
PrintOrDisp( "         and $numOfRv random Reverse sequences for each Forward sequence.\n");
PrintOrDisp("1st Selection\t: $numOfGoodCandidates was selected from $comb combination set.\n");
PrintOrDisp("2nd Selection\t: Fw$ii x Rv$jj was selected from 1st Selection set.\n");

################## Result of 1st selection
PrintOrDisp( "\n<Detail>\n" );
PrintOrDisp( "1st Selection (Score Board : Good candidates were showed in [ ])\n");
# label x
PrintOrDisp ("Fw\\Rv\t");
for(my $j=0; $j<$numOfRv; $j++) {
	my $jj = $j +1;
	PrintOrDisp ("%3d\t", $jj);
}
PrintOrDisp( "\n");

for(my $i=0; $i<$numOfFw; $i++) {
	#label y
	my $ii = $i +1;
	PrintOrDisp("%3d\t", $ii);
	for (my $j=0; $j<$numOfRv; $j++) {
		if ($Candidates[$i][$j] == 1) {
			PrintOrDisp("[%3d]\t", $firstScoreBoard[$i][$j]);
		}else{
			PrintOrDisp("%4d\t", $firstScoreBoard[$i][$j]);
		}
	}
	PrintOrDisp("\n");
}

################## show alignment of best primer pair
$flags{'visual'} = 'true';
my $annealedSeq = substr($primerFwSeq, 0, -$sizeOfOverLap);
$annealedSeq .= $primerRvSeq;
my $lenAnneal = length $annealedSeq;
my $lenFull = length $fullSequence;
my $lenFw = length $primerFwSeq;
my $lenRv = length $primerRvSeq;
my $offset = int($sizeOfOverLap/2);
my $partialScore;

# show full seq. alignment
PrintOrDisp("\nFw$ii ($lenFw bases) x Rv$jj ($lenRv bases)\n");
$score = getScore($primerFwSeq, $primerRvSeq); # show alignment this line
PrintOrDisp(" => Score $score\n");

# show partial seq. alignment
# partial Fw vs Rv
my $partialSeq = substr ($primerFwSeq, 0, -$offset);
$lenFw = length $partialSeq;
$lenRv = length $primerRvSeq;
PrintOrDisp("\nPartial_Fw$ii ($lenFw bases) x Rv$jj ($lenRv bases)\n");
$partialScore = getScore ($partialSeq, $primerRvSeq);
PrintOrDisp(" => Score $partialScore\n");
			
# Fw vs partial Rv
$partialSeq = substr ($primerRvSeq, $offset);
$lenFw = length $primerFwSeq;
$lenRv = length $partialSeq;
PrintOrDisp("\nFw$ii ($lenFw bases) x Partial_Rv$jj ($lenRv bases)\n");
$partialScore = getScore ($primerFwSeq, $partialSeq);
PrintOrDisp(" => Score $partialScore\n");

################## show sequence of best primer pair

# to complementary seq
$primerRvSeq =~ tr/ATGCatgc/TACGtacg/;
# to reverse seq
$primerRvSeq = reverse $primerRvSeq;

$lenFull = length $fullSequence;
$lenFw = length $primerFwSeq;
$lenRv = length $primerRvSeq;

PrintOrDisp("\nFull Sequence (length =$lenFull)\n$fullSequence\n\n");
PrintOrDisp("Suggested Forward Primer (length = $lenFw)\n\"$primerFwSeq\"\n\n");
PrintOrDisp("Suggested Reverse Primer (length = $lenRv)\n\"$primerRvSeq\"\n\n");
PrintOrDisp("Annealed Primers Sequence (length =$lenAnneal)\n\"$annealedSeq\"\n\n");
my $pepByReadingFrame0 = dna2peptide( $annealedSeq );
my $pepByReadingFrame1 = dna2peptide( substr($annealedSeq, 1) );
my $pepByReadingFrame2 = dna2peptide( substr($annealedSeq, 2) );
PrintOrDisp("Translations of Annealed Primers Sequence\n");
PrintOrDisp("Reading Frame0\n\"$pepByReadingFrame0\"\n");
PrintOrDisp("Reading Frame1\n\"$pepByReadingFrame1\"\n");
PrintOrDisp("Reading Frame2\n\"$pepByReadingFrame2\"\n");

################## For Debug ###
if($flags{'debug'} eq 'true') {
	PrintOrDisp("\nFor Debugging\n");
	PrintOrDisp("SortByScore1 \t\t\t\tSortByScore2\t\t\t\tTotalOrder\n");
	PrintOrDisp("Score\t\tID\tFw\tRv\tScore\t\tID\tFw\tRv\tOrder\tID\tTotalScore(2partial)\n");
	foreach (my $i=0; $i<@SortByScore1; $i++) {
		PrintOrDisp("$SortByScore1[$i]\t$SortByScore2[$i]\t$totalOrder[$i]\n");
	}
}


if ( $flags{'output'} eq 'true' ) {
	close OUTPUT_FILE;
}

##################################### End of Main


##################################### Sub functions
####################### SUB PrintOrDisp
sub PrintOrDisp {
	my (@args) = @_;
	my $format;
	my $string;
	my $isFormat = 'false';
	
	if (@args >1 ) {
		$format = shift;
		$string = shift;
		$isFormat = 'true';
	}else{
		$string = shift;
	}
	
	#print in file
	if ( $flags{'output'} eq 'true' ) {
		if ($isFormat eq 'true') {
			printf OUTPUT_FILE $format, $string;
		}else{
			printf OUTPUT_FILE $string;
		}
		
	#display in stdout
	}else{
		if ($isFormat eq 'true') {
			printf $format, $string;
		}else{
			printf $string;
		}
	}
}

####################### SUB getScore
sub getScore {

my (@arguments) = @_;
my $seq1 = $arguments[0];
my $seq2 = $arguments[1];

my ($selectingCandidate) = 'false';
my ($ctr1, $ctr2);
if( defined ($arguments[2]) && defined ($arguments[3]) ) {
	$selectingCandidate = 'true';
	$ctr1 = $arguments[2];
	$ctr2 = $arguments[3];
}

# scoring scheme
my $MATCH_AT    =  2;
my $MATCH_GC    =  3;
my $MISMATCH = -9;	#-9
my $GAP      = -13;	#-12
my $EXTGAP   = -2;

# initialization
my @matrix;
$matrix[0][0]{score}   = 0;
$matrix[0][0]{pointer} = "none";
for(my $j = 1; $j <= length($seq1); $j++) {
	$matrix[0][$j]{score}   = 0;
	$matrix[0][$j]{pointer} = "none";
}
for (my $i = 1; $i <= length($seq2); $i++) {
	$matrix[$i][0]{score}   = 0;
	$matrix[$i][0]{pointer} = "none";
}

# fill
my $max_i     = 0;
my $max_j     = 0;
my $max_score = 0;

for(my $i = 1; $i <= length($seq2); $i++) {
	for(my $j = 1; $j <= length($seq1); $j++) {
		my ($diagonal_score, $left_score, $up_score);
		
		# calculate match score
		my $letter1 = substr($seq1, $j-1, 1);
		my $letter2 = substr($seq2, $i-1, 1);		
		if ($letter1 eq $letter2) {
			if ($letter1 =~ /[AT]/i) {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH_AT;
			}else {
				$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH_GC;
			}
		}
		else {
			$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
		}
		
		# calculate gap scores
		if ( $matrix[$i-1][$j]{pointer} eq "up" ) {
			$up_score   = $matrix[$i-1][$j]{score} + $EXTGAP;
		}else {
			$up_score   = $matrix[$i-1][$j]{score} + $GAP;
		}
		if ( $matrix[$i][$j-1]{pointer} eq "pointer" ) {
			$left_score = $matrix[$i][$j-1]{score} + $EXTGAP;
		}else{
			$left_score = $matrix[$i][$j-1]{score} + $GAP;
		}
		
		if ($diagonal_score <= 0 and $up_score <= 0 and $left_score <= 0) {
			$matrix[$i][$j]{score}   = 0;
			$matrix[$i][$j]{pointer} = "none";
			next; # terminate this iteration of the loop
		}
		
		# choose best score
		if ($diagonal_score >= $up_score) {
			if ($diagonal_score >= $left_score) {
				$matrix[$i][$j]{score}   = $diagonal_score;
				$matrix[$i][$j]{pointer} = "diagonal";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		} else {
			if ($up_score >= $left_score) {
				$matrix[$i][$j]{score}   = $up_score;
				$matrix[$i][$j]{pointer} = "up";
			}
			else {
				$matrix[$i][$j]{score}   = $left_score;
				$matrix[$i][$j]{pointer} = "left";
			}
		}
		
		# set maximum score
		if ($matrix[$i][$j]{score} > $max_score) {
			$max_i     = $i;
			$max_j     = $j;
			$max_score = $matrix[$i][$j]{score};
		}
	}
}

my $align1 = "";
my $align2 = "";
my $matching = "";
my $no_matched = 0;
my $no_gap_sub = 0;
my $no_gap_que = 0;

my $j = $max_j;
my $i = $max_i;

while (1) {
	last if $matrix[$i][$j]{pointer} eq "none";
	
	if ($matrix[$i][$j]{pointer} eq "diagonal") {
		my $letter1 = substr($seq1, $j-1, 1);
		my $letter2 = substr($seq2, $i-1, 1);
		if ($letter1 eq $letter2) {
			if (substr($seq1, $j-2, 1) eq substr($seq2, $i-2, 1)
				&& substr($seq1, $j, 1) eq substr($seq2, $i, 1)) {
				$matching .= "|";
			}else{
				$matching .= ":";
			}
		}else{
			$matching .= " ";
		}
		
		$align1 .= $letter1;
		
		# change quiry to reverse-complementary seq
		$letter2 =~ tr/ACGTacgt/TGCAtgca/;
		$align2 .= $letter2;
		
		$no_matched++;
		$i--; $j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "left") {
		$align1 .= substr($seq1, $j-1, 1);
		$align2 .= "-";
		$matching .= " ";
		$no_gap_que++;
		$j--;
	}
	elsif ($matrix[$i][$j]{pointer} eq "up") {
		$align1 .= "-";
		
		# change quiry to reverse-complementary seq
		my $letter2 = substr($seq2, $i-1, 1);
		$letter2 =~ tr/ACGTacgt/TGCAtgca/;
		$align2 .= $letter2;
		
		$matching .= " ";
		$no_gap_sub++;
		$i--;
	}	
}

$align1 = reverse $align1;
$align2 = reverse $align2;
$matching = reverse $matching;

#$align1 = substr($seq1, 0, $j).$align1.substr($seq1, $max_j, length($seq1)-$max_j);
#$align2 = substr($seq2, 0, $i).$align2.substr($seq1, $max_i, length($seq2)-$max_i);

#print "$align1\n";
#print "$matching\n";
#print "$align2\n";

my $l = 0;
my @subject ='';
my @query ='';
my @matched = '';

while(1) {
    if(length $align1 < 50 * ($l+1)) {
	$subject[$l] = $align1;
	$matched[$l] = $matching;
	$query[$l] = $align2;
	#print "last\n";
	#print $subject[$l]."\n";
	last;
    }
    $subject[$l] = substr($align1, 0, 50);
    $matched[$l] = substr($matching, 0, 50);
    $query[$l] = substr($align2, 0, 50);
    $align1 =~ s/.{50}//;
    $matching =~ s/.{50}//;
    $align2 =~ s/.{50}//;
    $l++;
}

#statics
#my $identy = $no_matched / (50 * $l + length($subject[$l])) * 100;

#print result
my $line;
for($line=0; $line <= $l; $line++){
    $j = $j + 50 * $line;
    $i = $i + 50 * $line;
    
    if ($flags{"visual"} eq "true") {
        write;    
    }
    
    # added for check best matching in overlap
    # this should be removed if you want to use this sub-function at other programs
    if ($selectingCandidate eq 'true') {
        if ( $i==0 && length($query[$line]) == $sizeOfOverLap ) { # the start of reverse primer overlaps.
        	$Candidates[$ctr1][$ctr2] = 1; # yes
        }else{
        	$Candidates[$ctr1][$ctr2] = 0; # no
        }
    }
}

return $max_score;




format STDOUT =
Forward: @#### @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @####
$j, $subject[$line], length($subject[$line]) + $j
               @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
$matched[$line]
Reverse: @#### @<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< @####
$i, $query[$line], length($query[$line]) + $i
.

#format STDOUT_TOP =
#                     Allignment Reports
#  Match point: @##                       Identy: @#.## %
#$MATCH_AT, $identy
#     Mismatch: @##                  No. matched: @####
#$MISMATCH, $no_matched
#  Gap fanalty: @##             No. gap of sbjct: @####
#$GAP, $no_gap_sub
#      Ext-gap: @##             No. gap of query: @####
#$EXTGAP, $no_gap_que
#                                    Total score: @####
#$max_score
#Allignment:
#       Start                      Sequence                        End
#       --------------------------------------------------------------
#.
} # end of sub getScore()




####################### SUB rand_nucleotide
# Sub-function
# Random Nucleotide
# 2003, 11, 25 Kim Yu-sung. Exercise 7.7
#pass the probabilitis of each nucleotide into subroutine.
#and subroutines returns nucleotide.
# Warning: make sure you call srand to seed the
#	random number generator before you call this subroutine.
# usage: rand_nucleotide(pro1, pro2, pro3, pro4), pro1+pro2+pro3+pro4 == 1
sub rand_nucleotide {
	my(@probability) = @_;
    my($a, $c, $g, $t, $AT, $CG) = 0;
    my(@nucleotide) = qw/A C G T/;
    my($probability);

    if( !defined($probability[0]) ){
    	return $nucleotide[ rand(@nucleotide) ];
    }elsif( scalar(@probability) == 2 ){
		$AT = shift(@probability);
		$CG = shift(@probability);
		unless (($AT + $CG) == 1){
			die "not correct probability\n";
		}
	}elsif( scalar(@probability) == 4 ){
		($a,$c,$g,$t) = @probability;
		unless( ($a+$c+$g+$t) == 1 ){
			die "not correct probability\n";
		}
	}else{
		die "not correct probability\n";
	}

    $probability = rand(1);
    if( $AT ){
    	if($probability < $AT){
    		if( rand(1) < 0.5){
    			return 'A';
    		}else{ return 'T'; }
    	}else{
    		if( rand(1) < 0.5 ){
    			return 'C';
    		}else{ return 'G';}
    	}
    }else{
    	if($probability < $a){
    		return 'A';}
    	elsif($probability < $a+$c){
    		return 'C';}
    	elsif($probability < $a+$c+$g){
    		return 'G';}
    	else{
    		return 'T';
    	}
    }
}

####################### SUB rand_nucleotideXs
# These are the IUB ambiguity codes
# (Eur. J. Biochem. 150: 1-5, 1985):
# R = G or A
# Y = C or T
# M = A or C
# K = G or T
# S = G or C
# W = A or T
# B = not A (C or G or T)
# D = not C (A or G or T)
# H = not G (A or C or T)
# V = not T (A or C or G)
# N = A or C or G or T 

# N : A, C, G, T
sub rand_nucleotide_N {
	return rand_nucleotide (0.25, 0.25, 0.25, 0.25);
}
sub rand_nucleotide_weakN {
	return rand_nucleotide (0.4, 0.1, 0.1, 0.4);
}
sub rand_nucleotide_strongN {
	return rand_nucleotide (0.1, 0.4, 0.4, 0.1);
}

# Y : pYrimidine, C, T
sub rand_nucleotide_Y {
	return rand_nucleotide (0.0, 0.5, 0.0, 0.5);
}
sub rand_nucleotide_weakY {
	return rand_nucleotide (0.0, 0.2, 0.0, 0.8);
}
sub rand_nucleotide_strongY {
	return rand_nucleotide (0.0, 0.8, 0.0, 0.2);
}

# R : puRine, A, G
sub rand_nucleotide_R {
	return rand_nucleotide (0.5, 0.0, 0.5, 0.0);
}
sub rand_nucleotide_weakR {
	return rand_nucleotide (0.8, 0.0, 0.2, 0.0);
}
sub rand_nucleotide_strongR {
	return rand_nucleotide (0.2, 0.0, 0.8, 0.0);
}

# S : Strong, C, G
sub rand_nucleotide_S {
	return rand_nucleotide (0.0, 0.5, 0.5, 0.0);
}

# W : Week, A, T
sub rand_nucleotide_W {
	return rand_nucleotide (0.5, 0.0, 0.0, 0.5);
}

# H : A, C, T
sub rand_nucleotide_H {
	return rand_nucleotide (0.33333, 0.33333, 0.0, 0.33334);
}
# H : A, C, T
sub rand_nucleotide_weakH {
	return rand_nucleotide (0.4, 0.2, 0.0, 0.4);
}
# H : A, C, T
sub rand_nucleotide_strongH {
	return rand_nucleotide (0.1, 0.8, 0.0, 0.1);
}

# V : A, C, G
sub rand_nucleotide_V {
	return rand_nucleotide (0.33333, 0.0, 0.33333, 0.33334);
}
sub rand_nucleotide_weakV {
	return rand_nucleotide (0.4, 0.0, 0.2, 0.4);
}
sub rand_nucleotide_strongV {
	return rand_nucleotide (0.1, 0.0, 0.8, 0.1);
}

# D : A, G, T
sub rand_nucleotide_D {
	return rand_nucleotide (0.33333, 0.33333, 0.33334, 0.0);
}
sub rand_nucleotide_weakD {
	return rand_nucleotide (0.8, 0.1, 0.1, 0.0);
}
sub rand_nucleotide_strongD {
	return rand_nucleotide (0.2, 0.4, 0.4, 0.0);
}

# B : C, G, T
sub rand_nucleotide_B {
	return rand_nucleotide (0.0, 0.33333, 0.33333, 0.33334);
}
sub rand_nucleotide_weakB {
	return rand_nucleotide (0.0, 0.1, 0.1, 0.8);
}
sub rand_nucleotide_strongB {
	return rand_nucleotide (0.0, 0.4, 0.4, 0.2);
}

sub dna2peptide {

    my($dna) = @_;

    # Initialize variables
    my $protein = '';

    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= codon2aa( substr($dna,$i,3) );
    }

    return $protein;
}
#
# codon2aa
#
# A subroutine to translate a DNA 3-character codon to an amino acid
#   Version 3, using hash lookup

sub codon2aa {
    my($codon) = @_;

    $codon = uc $codon;
 
    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'L',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '_',    # Stop
    'TAG' => '_',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '_',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'V',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

    if(exists $genetic_code{$codon}) {
        return $genetic_code{$codon};
    }else{

            print STDERR "Bad codon \"$codon\"!!\n";
            exit;
    }
}
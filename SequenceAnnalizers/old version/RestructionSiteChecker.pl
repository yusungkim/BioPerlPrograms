####### SequenceGeneratorNRestructionSiteChecker.pl
####### 12.Jun.2013 YuSung Kim
####### usage: SequenceGeneratorNRestructionSiteChecker.pl [options] <Input_File> [-o Output_File] [-e Enzyme]
####### 하는일 : 프로그램안에 입력된 주어진 배열 안에서, 지정한 제한효소 인식 부위의 출현 갯수와 위치를 알려준다.
#######         배열에 멀티 코돈을 사용했을 때 유용하게 이용할 수 있도록 제한효소 인식부위의 출현 확률도 같이 알려준다.
####### 알고리즘 : 코드보면 쉽게 알 수 있으니 생략하겠다.

sub version {
    print "\n";
    print "#############################################################\n";
    print "##         Sequence Generator Version 1.2                  ##\n";
    print "##     .............................................       ##\n";
    print "##             All Copywrights are reserved by             ##\n";
    print "##           Yusung KIM, 2013.6.12 ~ 2013.06.13            ##\n";
    print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com    ##\n";
    print "#############################################################\n";
}
sub usage{
    print "\nusage: perl $0 [options] <Input File> [Output File]\n";
    print "ex)\n perl SequenceGeneratorNRestructionSiteChecker.pl -o result.txt -e KpnI -e HindIII\n";
    print "options:\n";
    print "	-v	-> display the process\n";
    print "	-d	-> display info for debug\n";
    print " -p  -> generate primers\n";
    print "	-V	-> display version\n";
    print "	-o	-> specify output file to save the results\n";
    print "		   ex) -o output.txt\n";
    print "	-e	-> specify restruction enzyme not to cut by\n";
    print "		   ex) -e KpnI -e BamHI\n";
    print "	-h/?	-> display this\n";
}


################################### MAIN #################################

#!/usr/bin/perl -w
use warnings;
use strict;

############################################## DATA INPUT #################
# The sequence which you want to create using PCR, the final product
#my ($PCRedSeq) = "NGCTCGACCATGGACGGGAATTCAGACVHKTATRNSCAGTGGTTGVNHVNHRRSGGACCGVNHTCTGGTCGTCCACCACCARNSGGARRSKCTRNSGGTRRTVHKTWTRNSVAWTGGMWGVNHVNHRRSGGTCCAVNHAGCGGACGTCCTCCTCCGTCTGGTRRSKCTRNSGGARRTVHKTATRNSCAATGGCTCVNHVNHRRSGGACCGVNHAGTGGTVRWCCTCCGCCARNSGGTRNSKCTRRSGGTRRTVHKTWYRNSVAWTGGTTAVNHVNHRRSGGACCGVNHTCAGGACGTCCACCTCCTRNSGGTRRSKCTRRSGGATCCATCGAAGCAGCACTCACTGCAG";
#my ($PCRedSeq) = "NGCTCGACCATGGACGGGAATTCAGACVHKTATRNSCAGTGGTTGVNHVNHRRSGGACCGVNHTCTGGTCGTCCACCACCARNSGGARRSKCTRNSGGTRRTVHKTWTRNSVAWTGGMWGVNHVNHRRSGGTCCAVNHTCCGGACGTCCTCCTCCGTCTGGTTCCGGARNSGGARRTVHKTATRNSCAATGGCTCVNHVNHRRSGGACCGVNHAGTGGTVRWCCTCCGCCARNSGGTRNSKCTRRSGGTRRTVHKTWYRNSVAWTGGTTAVNHVNHRRSGGACCGVNHAGCGGACGTCCACCTCCTRNSGGTRRSGGTGCTGCTGGATCCACGAGTCAGCAC";
my ($PCRedSeq) = "NCGAGCTCGACCATGGACGGGAATTCAGACVHKTATRNSCAGTGGTTGVNHVNHRRSGGACCGVNHTCTGGTCGTCCACCACCARNSGGARRSKCTRNSGGTRRTVHKTWTRNSVAWTGGMWGVNHVNHRRSGGTCCAVNHTCTGGACGTCCTCCTCCGTCTGGTTCCGGARNSGGARRTVHKTATRNSCAATGGCTCVNHVNHRRSGGACCGVNHAGTGGTVRWCCTCCGCCARNSGGTRNSKCTRRSGGTRRTVHKTWYRNSVAWTGGTTAVNHVNHRRSGGACCGVNHAGCGGACGTCCACCTCCTRNSGGTRRSGGTGCTGCTGGATCCACGAGTCAGCAC";


$PCRedSeq =~ s/\s//g;		#Remove spaces
$PCRedSeq =~ tr/a-z/A-Z/;	#Capitalize
$PCRedSeq =~ s/[^AGCTNYRMKSWHVDB]/ERROR!!NOT_AllowedBase/ig;	# Remove non-DNA code letter if any.

# The sequencing primers which anneal to the both end of PCRedSeq
my ($AmplifyPrimerOf5End) = "AGGAGGAATAATCCATGGAT GGG AAT TCA";
my ($AmplifyPrimerOf3End) = "GGA TCC GGTGGTTCATCTGGAG";
$AmplifyPrimerOf5End =~ s/\s//g;
$AmplifyPrimerOf5End =~ tr/a-z/A-Z/;
$AmplifyPrimerOf3End =~ s/\s//g;
$AmplifyPrimerOf3End =~ tr/a-z/A-Z/;

# Linker which is used for annealing of PrimerFw and PrimerRv
#my ($LinkerSeq) = "GTG AGC GTT GGT GCG GAA  GGC GTG";	# V-S-A-G-A-E-G-V
#my ($LinkerSeq) = "GTG GCG AGC ATT CTG GAA  GGC GTG";	# V-A-S-I-L-E-G-V
my ($LinkerSeq) = "AGT GGT CGC CCG CCA CCG";
$LinkerSeq =~ s/\s//g;		# Remove spaces
$LinkerSeq =~ tr/a-z/A-Z/;	# Capitalize ($LinkerSeq = uc $LinkerSeq;)

############################################## DATA END ###################


######## Parameters for input/output
my $input_file = '';
my $output_file = 'result.txt'; #default output file



######## Parameters for internal calculations
my ($sizeOfPCRedSeq) = length $PCRedSeq;

my ($sizeOf5End) = length $AmplifyPrimerOf5End;
my ($sizeOf3End) = length $AmplifyPrimerOf3End;

my ($sizeOfOverLap) = length $LinkerSeq;				# OverLap Size

my (%enzymeDataBase); #DB read from file, RestructionEnzymeList, ex) "KpnI" => "GGTACC"
my $numOfEnzymeDB;

my @enzymes = ("EcoRI", "BamHI");	# investigate these enzyme cut sites.

# 3차 배열, 각 요소는 지정된 배열의 지정된 효소에 의한 인식부위의 첫염기
my @cutPositions = ();		# 2-dimention	[Index_of_Enzymes][Index_of_Xth_Cut_of_Cuts]

# multibases
# elements which represent the probablilty of appearance
my @multibases = (
	#A, C, G, T
	[1, 0, 0, 0],	# A		index:0
	[0, 1, 0, 0],	# C
	[0, 0, 1, 0],	# G	
	[0, 0, 0, 1],	# T
	[0, 0, 0.5, 0.5],	# K		index:4
	[0.5, 0.5, 0, 0],	# M
	[0.5, 0, 0, 0.5],	# W
	[0, 0.5, 0.5, 0],	# S
	[0, 0.5, 0, 0.5],	# Y
	[0.5, 0, 0.5, 0],	# R
	[0, 0.3333, 0.3333, 0.3333],	# B index:10
	[0.3333, 0, 0.3333, 0.3333],	# D
	[0.3333, 0.3333, 0, 0.3333],	# H
	[0.3333, 0.3333, 0.3333, 0],	# V
	[0.25, 0.25, 0.25, 0.25]	# N		index: 14
);

my (%multibaseIndex) = (
	"A" => 0, "C" => 1,
	"G" => 2, "T" => 3,
	"K" => 4, "M" => 5,
	"W" => 6, "S" => 7,
	"Y" => 8, "R" => 9,
	"B" => 10, "D" => 11, 
	"H" => 12, "V" => 13,
	"N" => 14
);

#$multibases[ $multibaseIndex{'A'} ] = \@A;


######## temporalary data storages


####### CHECKING OPTIONS

#if( @ARGV == 0 ) {
#	version();
#	usage();
#	exit;
#}

# set flags as default
my (%flags) = (
"visual" => "false", ## donot display the process.
"debug"  => "false", ## donot display info for debug.
"output" => "false",
"make primer"	=> "false"
);

if( defined($ARGV[0]) && ($ARGV[0] =~ /^-\S*V\S*/) ) {
    version();
    exit;
}
if ( defined($ARGV[0]) && ($ARGV[0] =~ /^-\S*[h?]\S*/) ) {
	usage();
	exit;
}

my $temp = '';
my $enzymeSpecifiedAlready = 'false';
while( defined($temp = shift) ) {
    
    # -j : just test given random sequence in this file
	if($temp =~ /^-\S*p\S*/){
		$flags{"make primer"} = "true";
	}
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
		my $enzymeName = shift;
		if($enzymeSpecifiedAlready eq 'false') {
			@enzymes = ();
			$enzymeSpecifiedAlready = 'true';
		}
		# check if the enzname fllows.
		if ((!defined $enzymeName) || ($enzymeName =~ /^-/)) {
			print "Wrong option usage:\n";
			print "\tEnzyme name should be followed after -e option\n";
			exit;
		}
		push @enzymes, $enzymeName;
	}
	# Wrong options
	if($temp =~ /^-[vdope]*([^dvope\d])+[dvpoe]*$/){
		print "Wrong option used: $1\n";
		usage();
		exit;
	}
}

################################################### PROCESS START #############
print "\n";

####################### read restruction enzyme from file
PrintOrDisp ("Processing .... ");
open(RES, "RestructionEnzymeList.txt") || die "Cannot open:\n\t$!";
while ( <RES> ) {
    chomp;
    $_ =~ /^(\S+)\s*,\s*(\S+)\s*/;
    my $name = $1;
    my $seq = $2;
    $seq =~ s/[^a-zA-Z]//;
    $enzymeDataBase{$name} = $seq;
    #print "$name => $enzymeDataBase{$name}\t";
}
close(RES);
$numOfEnzymeDB = keys %enzymeDataBase;
PrintOrDisp ( "$numOfEnzymeDB Restruction Enzymes were read.\n");


####################### pre-process
# remove everything excetp DNA
#$fullSequence =~ s/[^AGCTNYRSWHVDB]//ig;

# remove spaces
#$fullSequence =~ s/[\s]//g;

# remove numbers
#$fullSequence =~ s/\d//g;

# remove non-alphabet
#$fullSequence =~ s/[^A-Za-z]//g;

# change to Upper case
##$fullSequence = uc $fullSequence; # $fullSequence =~ tr/[a-z]/[A-Z]/;



####################### Restruction Enzyme Site Check
# 작성한 랜덤소팅된 멀티코돈 배열에서 지정한 제한효소 부위가 등장하는지 검사한다.
# 부위의 갯수가 없는 배열조합을 선택출력한다.
PrintOrDisp ("Processing .... ");


if ( $flags{'visual'} eq 'true' ) {
	my $seqByThreeLetters = GetSeqByThreeLetters($PCRedSeq);
	PrintOrDisp ("$seqByThreeLetters\n");
}
	
# Find restruction enzyme recognition site
@cutPositions = ();
my $cutNumberOfAllEnzyme = 0;
my $pOfNoCutAgainstAllEnzyme = 1.0;	# init for calc.
my $numOfCutsAgainstAllEnzyme = 0;	# init for calc.

# check for all specified restruction enzymes
for(my $j=0; $j<@enzymes; $j++) {
	my ($enzyme) = $enzymes[$j];
	my ($recognitionSite) = $enzymeDataBase{$enzyme};
	my $sizeOfRecognitionSite = length $recognitionSite;
	
	my $probability = 1.0;			# p of Cut
	my $pOfNoCutAtSinglePosition = 1.0;
	my $pOfNoCut = 1.0;			# init for calc.
	my $numOfCuts = 0;
	
	if ( $flags{'visual'} eq 'true' ) { PrintOrDisp ("$enzyme: \t$recognitionSite\n"); }
	
	# KMP 방법을 쓰면 더 빠를 수도 있겠지만, base가 4개밖에 없는 DNA배열특성상 효율이 크게 향상되지는 않을 것 같다.
	# 그냥 알기 쉽게 차근차근 비교하는 방법으로 search하자.
	# 먼저 DNA sequence를 array에 넣어 인덱스 조작을 쉽게 한도록한다.
	my (@PCRedSeqInArray) = split (//, $PCRedSeq);
	my (@recognitionSiteInArray) = split (//, $recognitionSite);
	
	# Find restruction enzyme recognition site
	# 배열의 처음부터 끝까지 비교
	for (my $ctr1 = 0; $ctr1 < $sizeOfPCRedSeq - $sizeOfRecognitionSite; $ctr1++) {
	
		# 지정된 효소의 인식부위와 일대일비교
		$probability = 1.0;
		for (my $ctr2 = 0; $ctr2 < $sizeOfRecognitionSite; $ctr2++) {
			my $DNA1 = $PCRedSeqInArray[$ctr1 + $ctr2];
			my $DNA2 = $recognitionSiteInArray[$ctr2];
			
			my $p = GetPOfSameBaseOccurance($DNA1, $DNA2); # single base
			
			if ( $flags{'debug'} eq 'true' ) { print "$p\t"; } # for debugging
			
			$probability *= $p;	# get probility of enzyme cut occurance at one position
			
			if($probability == 0) {
			# No Recognition in this region.
				$ctr2 = $sizeOfRecognitionSite;	# End the roop
			
				if ( $flags{'debug'} eq 'true' ) {
					print "\n";
				} # for debugging
			
			}else{
				# Do nothing just Check the next base
			}
		}
			
		if($probability != 0){ # if recognition site was found)
	
			# record cut positions at two-dim-array[Index_Enz][Index_Cut]
			$cutPositions[$j][$numOfCuts] = $ctr1;
			
			# increase cut count
			$numOfCuts++;
			
			# cacl. probability that enzyme cut never occurs
			# only consider the region between two amplify primers.
			$pOfNoCutAtSinglePosition = (1.0 - $probability);
			if($sizeOf5End < $ctr1 + $sizeOfRecognitionSite && $ctr1 < $sizeOfPCRedSeq - $sizeOf3End) {
				$pOfNoCut *= $pOfNoCutAtSinglePosition;
			}
			 
			if ( $flags{'visual'} eq 'true' ) {
				my $seq = substr ( $PCRedSeq, $ctr1, $sizeOfRecognitionSite );
				PrintOrDisp( "\tcut at %03d\t", $ctr1);
				PrintOrDisp( "\"$seq\"\t");
				PrintOrDisp( "p=%.6f\t", $probability);
				PrintOrDisp( "pOfNoCut=%.4f\n", $pOfNoCut);
			} # end of if
			
		} # end of if
		
	} # end of for() # finding cut site


	my $noCutPercentage = $pOfNoCut * 100;
	if ( $flags{'visual'} eq 'true' ) {
		PrintOrDisp (" => Possible cuts of %8s", $enzymes[$j]);
		PrintOrDisp (": $numOfCuts\tpNoCut: $noCutPercentage %%\n");
	}

	$numOfCutsAgainstAllEnzyme += $numOfCuts;
	$pOfNoCutAgainstAllEnzyme *= $pOfNoCut;

} # end of searching recognition site for each Enzyme
	
#@cutPositions = sort numSort @cutPositions;
	
# for debugging
if ( $flags{'debug'} eq 'true') {
	&read_array(\@cutPositions); # show the whole array
}


my $dummy_howMany = scalar @enzymes;
PrintOrDisp ("$dummy_howMany enzymes (");
for(my $i=0; $i<$dummy_howMany-1; $i++) { PrintOrDisp ("$enzymes[$i], "); }
PrintOrDisp ("$enzymes[$dummy_howMany-1])were tested.\n");

######################## Sort for Selection of best sequence
# sort seqs by High pOfNoCut
# sort algorithm : insertion sort

# calc. num of Frame shift cuts

my $numOfFrameShiftCuts = 0;
for (my $j = 0; $j < @enzymes; $j++) {
	for (my $k = 0; $k < $#{$cutPositions[$j]} + 1; $k++) {
		if ($cutPositions[$j][$k] % 3 != 0) {
			$numOfFrameShiftCuts++;
		}
	}
}


#################################### Report #############################################
PrintOrDisp( "\n<Report>\n" );
PrintOrDisp( "=============================================================================\n" );
#PrintOrDisp("\n");

my $seqByThreeLetters = '';
my @cutScreens = ();
	
# print cut statistics
my $numOfEnzymes = $#enzymes + 1;	#get size of an array, or you can do 'scalar @enzymes'
# $numOfEnzymes = @enzmyes;
# $numOfEnzymes = scalar @enzymes;
	
for(my $j=0; $j<@enzymes; $j++) {
	my $enzymeNum = $j+1;
	
	my $cuts = $#{$cutPositions[$j]} + 1;
	#my $cuts = $#{$cutPositions[$i]} + 1;
	
	# ex) #1 NdeI (6), #2 EcoRI (5), #3 ...
	PrintOrDisp ("#$enzymeNum:$enzymes[$j]($cuts)");
	
	if($j<$numOfEnzymes-1) {
		PrintOrDisp (" + ");
	}
}
# print total cuts
PrintOrDisp ("\t Total %d Cuts", $numOfCutsAgainstAllEnzyme);

# print total cuts which induce frame shift
PrintOrDisp ("(%d Frame Shift)", $numOfFrameShiftCuts);

# print posibillity of no-cut in the region between two amplfy primers.
PrintOrDisp ("\tRandom Region Cut = %.2f%%\n", (1.0-$pOfNoCutAgainstAllEnzyme)*100);

# print sequences
PrintOrDisp ("\n  [DNA Sequence] $sizeOfPCRedSeq base\n");
$seqByThreeLetters = GetSeqByThreeLetters($PCRedSeq);
PrintOrDisp ("          $seqByThreeLetters\n\n");

for (my $j = 0; $j < @enzymes; $j++) {
	my $cutScreen = '';		# cut position will be set on this string
	my $nextPosition;
	my $indexNext = 0;		# index for @totalCutPositions

	# set first cutPosition if there is.
	if( defined $cutPositions[$j][0]) {
		$nextPosition = $cutPositions[$j][0];

		# print every cut position on the top of seq
		for (my $k = 0; $k< $sizeOfPCRedSeq; $k++) {
			if($k<$nextPosition) {
			$cutScreen .= "-";
			}elsif($k == $nextPosition){
				my $EnzymeSite = $enzymeDataBase{$enzymes[$j]};
				
				if(++$indexNext < $#{$cutPositions[$j]}+1) {
					$nextPosition = $cutPositions[$j][$indexNext];
					if($k <= $nextPosition - length $EnzymeSite ) {
						$cutScreen .= $EnzymeSite;
						$k = $k - 1 + length $EnzymeSite;
					}else{
						$cutScreen .= "|";	# if another site exist closly behind
					}
				}else{
					$cutScreen .= $EnzymeSite;
					$nextPosition = $sizeOfPCRedSeq;	# no more cut
					$k = $k - 1 + length $EnzymeSite;
				}
			}
		}
		
		# print cut position by three letter spaced way
		$seqByThreeLetters = GetSeqByThreeLetters($cutScreen);
		my $enzymeNameForPrint = "[" . $enzymes[$j] . "]";
		PrintOrDisp ("%9s ", $enzymeNameForPrint);
		PrintOrDisp ("$seqByThreeLetters\n\n");

	# if theris no cut.
	}else {
		print"No cuts!\n";
	}
}	

# make primers
if ( $flags{'make primer'} eq 'true') {
	my $primerFw;
	my $primerRv;
	my $sizeOfFw;
	my $sizeOfRv;
	
	if( $PCRedSeq =~ m/$LinkerSeq/ ) {
	 	$primerFw = $` . $&;	# 5' + Annealing region
		$primerRv = $& . $';	# Annealing region + 3'
		$sizeOfFw = length $primerFw;
		$sizeOfRv = length $primerRv;
	}else{
		print "		Fatal Error	: Cannot find Annealing seq.\n";
	}
	
	# get reverse-complement seq
	$primerRv = GetReverseComplement($primerRv);
	
	PrintOrDisp ( "[Primer Fw] %d bases\n", $sizeOfFw );
	PrintOrDisp ( "%s\n", GetSeqByThreeLetters($primerFw) );
	PrintOrDisp ( "[Primer Rv] %d bases\n", $sizeOfRv);
	PrintOrDisp ( "%s\n\n", GetSeqByThreeLetters($primerRv) );
}

if ( $flags{'output'} eq 'true' ) {
	close OUTPUT_FILE;
}

##################################### End of Main ##############################

##################################### Sub functions ############################
####################### SUB GetReverseComplement
sub GetReverseComplement {
	my (@args) = @_;
	my $original = shift;
	
	my %complement = (
		  'G'=>'C','A'=>'T','T'=>'A','C'=>'G',
		  'N'	=>	'N',
		  
		  'W'	=>	'W',
		  'S'	=>	'S',
		  
		  'Y'	=>	'R',
		  'R'	=>	'Y',
		  'K'	=>	'M',
		  'M'	=>	'K',
		  
		  'V'	=>	'B',
		  'B'	=>	'V',
		  'H'	=>	'D',
		  'D'	=>	'H'	);
	
	## reverse
	my $reversed = reverse $original;
	
	## complement
	my (@dummy) = split (//, $reversed);
	for(my $i=0; $i<@dummy; $i++) {
		$dummy[$i] = $complement{$dummy[$i]};
	}
	
	return join ('', @dummy);
}

####################### SUB GetSeqByThreeLetters
sub GetSeqByThreeLetters {
	my (@args) = @_;
	my $seq = shift;
	my $length = length $seq;
	
	if ($length < 1) {
		return;
	}
	
	my $seqByThreeLetters = substr($seq, 0, 1);
	
	for(my $i=1; $i<$length; $i++) {
		if ($i % 3 == 0) {
			$seqByThreeLetters .= ' ';
		}
		$seqByThreeLetters .= substr($seq, $i, 1);
	}
	
	return $seqByThreeLetters;
}

####################### SUB GetPOfSameBaseOccurance
sub GetPOfSameBaseOccurance {
	my (@args) = @_;
	my $multibase = shift;
	my $base = shift;

	my $indexMultibase = $multibaseIndex{ $multibase };
	my $indexBase = $multibaseIndex{ $base };

	my $probability = $multibases[$indexMultibase][$indexBase];
	
	if ( $flags{'debug'} eq 'true' ) {
		print "\'$multibase\'";
	}
	
	return $probability;
}

####################### SUB TestTranslation
sub TestTranslation {
	my (@args) = @_;
	my $multicodonSeq;
	my $repeat = 5;			# Number of sampling from given multi-codon sequence. Default = 2
	
	# copy PCRedSeq
	if (@args > 1) {
		$multicodonSeq = shift;
		$repeat = shift;
	}else{
		$multicodonSeq = shift;
	}
	
	#PrintOrDisp("Translations of Annealed Primers\n");
	#PrintOrDisp("Multicodon DNA Sequence:\n\t$multicodonSeq\n");
	
	for(my $ctr=0; $ctr < $repeat; $ctr++) {
		my $sampleDNASeq = $multicodonSeq;
	
		$sampleDNASeq =~ s/R/&rand_nucleotide_R/ge; # puRine		: A, G
		$sampleDNASeq =~ s/Y/&rand_nucleotide_Y/ge; # pYrimidine	: C, T
		$sampleDNASeq =~ s/S/&rand_nucleotide_S/ge;	# Strong		: C, G
		$sampleDNASeq =~ s/W/&rand_nucleotide_W/ge; # Weak			: A, T
		$sampleDNASeq =~ s/M/&rand_nucleotide_M/ge; # aMino			: A, C
		$sampleDNASeq =~ s/K/&rand_nucleotide_K/ge; # Keto			: G, T
		$sampleDNASeq =~ s/H/&rand_nucleotide_H/ge; # H				: T, C, A
		$sampleDNASeq =~ s/V/&rand_nucleotide_V/ge; # V				: A, C, G
		$sampleDNASeq =~ s/D/&rand_nucleotide_D/ge; # D				: A, G, T
		$sampleDNASeq =~ s/B/&rand_nucleotide_B/ge; # B				: C, G, T
		$sampleDNASeq =~ s/N/&rand_nucleotide_N/ge; # N				: A, C, G, T
	
		my $pepByReadingFrame0 = dna2peptide( $sampleDNASeq );
	
		#PrintOrDisp(" $ctr) $sampleDNASeq\n");
		PrintOrDisp(">Pep%d\n", $ctr+1);
		PrintOrDisp("         $pepByReadingFrame0\n");
	}
}

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

# M : aMino, A or C
sub rand_nucleotide_M {
	return rand_nucleotide (0.5, 0.5, 0.0, 0.0);
}

# K : Keto, G, T
sub rand_nucleotide_K {
	return rand_nucleotide (0.0, 0.0, 0.5, 0.5);
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
	return rand_nucleotide (0.33333, 0.33333, 0.33334, 0);
}
sub rand_nucleotide_weakV {
	return rand_nucleotide (0.8, 0.1, 0.1, 0.0);
}
sub rand_nucleotide_strongV {
	return rand_nucleotide (0.2, 0.4, 0.4, 0.0);
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

sub numSort {
  if ($a < $b) { return -1; }
  elsif ($a == $b) { return 0;}
  elsif ($a > $b) { return 1; }
}
sub read_array {
    my $array_reference = shift;

    foreach (@$array_reference) {
        if (ref($_) eq 'ARRAY') {
            &read_array($_);
            print "\n";
        } else {
            print $_ . " ";
        }
    }
}
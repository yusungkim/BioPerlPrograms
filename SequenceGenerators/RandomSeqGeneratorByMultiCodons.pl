####### RandomSeqGeneratorByMultiCodons.pl
#######  upgraded from "SequenceGeneratorNRestructionSiteChecker.pl"
####### 7.Apr.2012 YuSung Kim
####### usage: RandomSeqGeneratorByMultiCodons.pl [options] [-o Output_File] [-e Enzyme]
####### 하는일 : 입력된 Multi-codons과 사용횟수를 이용하여 Multi-codons을 랜덤으로 배치시킨 배열을 얻어낸다.
#######         업그레이드된 내용
#######         1. 두개의 긴 랜덤 프라이머를 anneal하여 랜덤시퀀스 이중나선을 만드는 방법을 포기한다.
#######         대신 길이가 다른 (약 두배) 랜덤 프라이머 두개를 각각 PCR을 통하여 증폭, 2중나선화 시키고,
#######         제한효소 하나를 사용하여 PCR산물을 연결시킨다. 길이가 다른 PCR산물이기에 hetero dimer, 즉,
#######         원하는 목적 랜덤시퀀스 이중나선을 만들 수 있다.
#######         이 방법이 anneal을 통하여 제작하는 방법보다 낳은 점은 specific한 anneal의 어려움으로 인한 부산물의
#######         생성을 낮출 수 있다 (2중나선화 작업때에 사용되는 길이가 더 작기때문에)
#######         2. 위와 같은 내용을 실현하기 위해, 주어진 멀티시퀀스와 그 사용빈도를 사용하여 랜덤프라이머를 산출한다.
#######         3. 싱글밴드를 만들기위한 성질이 가장높은 프라이얼러티 (예상anneal을 산출해서 제일 낮은 것을 선택)

####### 제작상의 아이디어들
####### 1. PCR 산물의 결합에 사용할 restriction enzyme은 SalI을 사용.
#######     재고가 있고, 높은 농도이며 (15U/ul), 선택한 멀티코돈조합에서 잘 발생하지 않는 인식부위라서(XhoI도 그러함)
####### 2. 

####### 알고리즘 : 코드보면 쉽게 알 수 있으니 생략하겠다.

sub version {
    print "\n";
    print "#############################################################\n";
    print "##         Sequence Generator Version 2.0                  ##\n";
    print "##     .............................................       ##\n";
    print "##             All Copywrights are reserved by             ##\n";
    print "##           Yusung KIM, 2012.2.23 ~ 2012.04.6             ##\n";
    print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com    ##\n";
    print "#############################################################\n";
}
sub usage{
    print "\nusage: perl $0 [options] [-o Output_File] [-e Enzyme]\n";
    print "ex)\n perl $0 -o result.txt -e KpnI -e HindIII\n";
    print "options:\n";
    print "	-v	-> display the process\n";
    print "	-d	-> display info for debug\n";
    print "	-V	-> display version\n";
    print "	-p	-> display primer pairs\n";
    print "	-j	-> Just test given full seq in this file\n";
    print "	-s	-> Get Sequences as specified number\n";
	print "        ex) -s8\n";
	print "	-t	-> display possible translation\n";
	print "        ex) -t5\n";
    print "	-o	-> specify output file to save the results\n";
    print "		   ex) -o output.txt\n";
    print "	-e	-> specify restruction enzyme not to cut by\n";
    print "		   ex) -e KpnI -e BamHI\n";
    print "	-a	-> codon usage anysis from given multi-codon set\n";
    print "	-h/?	-> display this\n";
}


################################### MAIN #################################

#!/usr/bin/perl -w
use warnings;
use strict;

# Specify the best random-sequence for primer making

# Shield_F1 #1
#my $SequenceForjustTestGivenFullSeq = "GAW TKK RNS TKK RNS VAW TKK RMA GAA RNS RRT GAA GAW SMR VHG RWG SMR RNS SMW VHG RMA DHT VMW RRS TSS GAW RWG RNS GAA RNS DHT VAW VHG TSS KRT GND RRS SWG VMW SMW GAW WNY CTC TCA GCC TAC GAG GCC TTA CCC TTC VMW SMW DHT GAW RWG KRT RRT KRT GAA GAA RRS TSS SWG GAW RMA RRS SMR WNY GAA VHG GND RRT RMA RNS GND";
#my $SequenceForjustTestGivenFullSeq ="GAW TKK RNS TKK RNS VAW TKK RMA GAA RNS RRT GAA GAW SMR VHG RWG SMR TGT RNS SMW VHG RMA GAA DHT VMW RRS TSS GAW RWG RNS GAA RNS DHT GAA VHG TSS KRT GND RRS GAA SWG VMW SMW TGC GAW WNY CTC TCA GCC TAC GAG GCC TTA CCC TTC VMW SMW DHT GAW RWG KRT RRT KRT GAA GAA RRS TSS SWG GAW RMA RRS GAA SMR WNY GAA TGT VHG GND RRT RMA RNS TGT GND";
#my $SequenceForjustTestGivenFullSeq = "GAW TKK RNS TKK RNS VAW TKK RMA GAA RNS RRT GAA GAW SMR VHG RWG SMR TGT RNS TGT VHG RMA GAA DHT VMW RRS TSS GAW RWG RNS GAA RNS DHT GAA VHG TSS KRT GND RRS GAA SWG RAG SMW TGC GAW WNY CTC TCA GCC TAC GAG GCC TTA CCC TTC VMW SMW DHT AAT GAW RWG KRT RRT KRT GAA GAA RRS TSS SWG GAW RMA RRS GAA GAC SMR WNY GAA TGT VHG GND RRT TTT RMA RNS TGT GND";
#my $SequenceForjustTestGivenFullSeq = "GAW TKK RNS TKK RNS VAW TKK RMA GAA RNS RRT GAA GAW SMR VHG RWG SMR TGT RNS TGT VHG RMA GAA DHT VMW RRS TSS GAW RWG RNS GAA RNS DHT GAA VHG TSS KRT GND RRS GAA SWG RAG SMW TGC GAW WNY CTC TCA GCC TAC GAG GCC TTA CCC TTC VMW SMW DHT AAT GAW RWG KRT RRT KRT KGT GAA RRS TSS SWG GAW RMA RRS GAA GAC SMR WNY GAA VHG TGT GND RRT YTC RMA GND TGT RNS";

# Shield_F1 #1 x #2
#my $ = "NNY GND RMA VNS GND VHG RAD GND VHG NBB SWG RMA VNS VHG VNS NBB DHT SWG SWG SWG SWG VHG RMA NNY GNR RAD SWG NBB NNY VNS NNK DHT RAD VHG SWG VHG GNR GNR RMA GNY VNS NBB SWG SRY GNR GTG AGC GTT GGT GCG GAA GGC GTG NNY GNR VNS GND DHT NNY RMA GNY DHT GNY GND VHG DHT SWG VHG VHG GND VNS VNS VNS DHT VNS SWG VNS DHT DHT RAD NNK VNS SWG NBB SWG RAD VHG DHT DHT GNR NNY RAD WNY VNS NNY SWG NNY RMA";
#my $SequenceForjustTestGivenFullSeq = "RMA GAW VMW RMA TKK GAC KRT TSS DHT RNS DHT RRT VHG DHT DHT RMA GAW SMW TGT VHG RNS TSS DHT RNS WNY RRS TGT TGC SMW SMR RWG SWG GAC VHG DHT GAW TGC RRT VMW CTC TCA GCC TAC GAG GCC TTA CCC TTC TKK TGC RNS GAW GND VAW SMR GAA KRT SMR RMA VAW RMA SWG GND RRS KRT GAA RWG GAW RRS WNY RWG VMW GND RNS TSS GAA SMW GAA GAW RRT VHG RNS TKK VMW RRS RNS RRS";
#my $SequenceForjustTestGivenFullSeq = "AAT ACC ATG GAA CAT ATG VNS RRS VMW VHK DHT RNS RNS DHT RNS RNS DHT VHK VHK SNN VHG RMA NNT SMR SMN SWG VMW RMA VHK VMW RMA RRS VMW RRS SNN VNS GND SNN RRS NBB VNS RRS VHK WNY RNS VMW GND RMA RNS NNT SVN VMW SMR NNT NNT RNS VHG SMN RMA RRS VAW RMA RNS NNT GTG AGC GTT GGT GCG GAA GGC GTG SWG RRS DHT WNY RMA VHK GND SMR RMA RNS NNT GND RRS SVN VNS NNT RRS NNT SMN WNY VNS DBK RNS NNT VHK NBB SVN VHG VAW VHK SMR SMN VAW NNT SNN VHK NBB NNT SVN VHG SWG GND DBK SWG VHK SMR SMN RNS DHT VHG RMA RNS DBK NNT SNN RRS SVN RMA GAA TTC AAG CTT ATA";
#my $SequenceForjustTestGivenFullSeq = "RRS NNT RMA VHK RNS NNT RRS VHK SWG VHK SMN DBK WNY RRS RRS VHK SWG DHT RMA SVN NNT VNS SMN RRS NBB DHT VNS SMN NNT SMN SMR VHG VMW GND RMA NBB RMA RMA SMR RRS VHG VMW RNS NNT RNS DHT SNN VMW VHG DBK VHK SVN SMR RMA RNS RNS GTG AGC GTT GGT GCG GAA GGC GTG GND DBK NNT RNS RRS GND VHK DHT VMW SVN RMA VHG GND RNS RMA SNN SMR SWG VHK RNS NNT SVN WNY RRS SMR VAW NNT DHT NNT RMA SMN GND NNT NBB NNT VHK SNN VHK VNS RNS SNN VMW VHK SWG SVN VAW RRS SNN RNS VAW RRS VHG RMA VNS NNT VNS";
#my $SequenceForjustTestGivenFullSeq = "DHT GND RMA VAW RRS SWG SVN SMR VMW VNS NNT RMA GND VMW VMW SVN SWG NNT VNS SNN SMN VHK VHK RRS SMN VHK SNN VAW RRS SNN GND SMR SVN VHG VAW NNT RNS SNN RRS NNK NNT SNN RRS RMA RNS RMA RMA VMW SMR SMN SVN SMN RNS NNT RRS SMN SWG VAW GTG AGC GTT GGT GCG GAA GGC GTG VMW NNT VHK RMA RRS VHK VNS VHK DHT VNS DHT NNT NNT SWG RMA SMR VHG SWG RNS SWG VNS RNS VHK RRS RNS RRS NNT NNT DHT SNN RRS VAW RMA GND SVN VHK GND RMA VHK DHT VHG RRS SVN RMA RMA VMW VNS VHG VAW SMR NNT VHG SMR VHK SMN VHK VHG GND DHT";
#my SequenceForjustTestGivenFullSeq = "VHG RMA RMA VAW RMA VMW VHK DHT SMN VNS NNT NNT NNT VAW RRS VHG RNS DHT NNT RRS SNN VMW VHG VHK SMR RMA VAW SWG SMN NNK RRS SNN VHK VHK VNS GND SNN VMW SNN SVN RNS RMA VHK DHT SNN RNS SWG DHT SVN NNT DHT NNT SMR SMN RNS NNT GTG AGC GTT GGT GCG GAA GGC GTG VMW GND VAW DHT VMW GND NNT SMR SVN VHK SMR VAW VMW RRS VHG GND SVN NNT VHK VNS NNT SWG RMA VNS RMA SMN VHK SVN RRS RMA VNS SMN GND SWG VHG VHK RNS RMA SNN SMR SMN VHK RMA NNT SMR SWG RNS VHG GND VNS SVN SWG RMA VAW VHK RRS";
#my $SequenceForjustTestGivenFullSeq = "AAT ACC ATG GAA CAT ATG NNY GND RMA VNS GND VHG RAD GND VHG NBB SWG RMA VNS VHG VNS NBB DHT SWG SWG SWG SWG VHG RMA NNY GNR RAD SWG NBB NNY VNS NNK DHT RAD VHG SWG VHG GNR GNR RMA GNY VNS NBB SWG SRY GNR GTG AGC GTT GGT GCG GAA GGC GTG NNY GNR VNS GND DHT NNY RMA GNY DHT GNY GND VHG DHT SWG VHG VHG GND VNS VNS VNS DHT VNS SWG VNS DHT DHT RAD NNK VNS SWG NBB SWG RAD VHG DHT DHT GNR NNY RAD WNY VNS NNY SWG NNY RMA GAA TTC AAG CTT ATA";

my $SequenceForjustTestGivenFullSeq = "SMW VYG RRS WNY RRS GAW RWG VYG VHG SMR RWG MDT RRT VMW GAA DHT VHG VHG MDT SMW RRT VHG VMW TGT DHT RMA RRT VMW RMA SMW SMW TKK RMA SWG RRS RMA NYK RRS DHT KRT TSS GND TKK SWG RMA GAW TSS CTG TCT AAA GAC CCT AAT GAG AAA CGC GAC RWG WNY KRT MDT GAW DHT DHT DHT RMA MDT VYG GAA RRS VMW VAW GND GAA GAT RNS GND GAW RRS RMA RRS VHG SMR RRS RMA GAW DHT SMW GAW TKK VAW CAG RMA MDT VHG TSS GAA RRS DHT SMR GAA ATT DHT KRT";


# Linker which is used for annealing of PrimerFw and PrimerRv
#my $LinkerSeq = "CTC TCA GCC TAC GAG GCC TTA CCC TTC";	# L-S-A-Y-E-A-L-P-F
my $LinkerSeq = "CTG TCT AAA GAC CCT AAT GAG AAA CGC GAC"; #L-S-K-D-P-N-E-K-R-D, GFP loop
#my ($LinkerSeq) = "GTG AGC GTT GGT GCG GAA GGC GTG";	# V-S-A-G-A-E-G-V
#my ($LinkerSeq) = "GTG GCG AGC ATT CTG GAA GGC GTG";	# V-A-S-I-L-E-G-V

$LinkerSeq =~ s/\s//g;	# Remove spaces
$LinkerSeq =~ tr/[a-z]/[A-Z]/;	# To Upper Character

my ($sizeOfOverLap) = length $LinkerSeq;	# OverLap Size

# primers
my @primersOfFw = ();
my @primersOfRv = ();

# default cloning sites
my $cloningSite5 = "tCC ATG GAA TGG ATC CAA GAA GGA GAT";	# 5'-NcoI-BamHI-3'
my $cloningSite3 = "GAA CTT GTA TCT GTC GAA TTC GAA AA";		# 5'-EcoRI-3'
#my $cloningSite5 = "CCATGGAACATATG";	# 5' - NcoI - AA - NdeI - 3', frame of NcoI is diff with NdeI
#my $cloningSite3 = "GAATTCAAGCTT";		# 5' - EcoRI - HindIII - 3'

$cloningSite5 =~ s/\s//g;	# Remove spaces
$cloningSite3 =~ s/\s//g;	# Remove spaces
$cloningSite5 =~ tr/[a-z]/[A-Z]/;	# to upper character
$cloningSite3 =~ tr/[a-z]/[A-Z]/;	# to upper character

# possible annealing positions
my %annealSeqs;				# hash {'forward'/'reverse'} = sequence
my @annealPositions = ();	# 3-dimension	[Index_Seq][Index_Enz][Index_Anneal]
my @annealLengthes = ();	# 3-dimension	[Index_Seq][Index_Enz][Index_Anneal]
my @annealPossibility = ();	# 3-dimension	[Index_Seq][Index_Enz][Index_Anneal]
my @pOfNoMatchAllSeqs = ();	# 1-dimension	[Index_Seq]
my @numOfMatchesAllSeqs=();	# 1-dimension	[Index_Seq]

my (%multicodons) = (
"DBK"  => 0,
"DHT"   => 9,
"GAW"	=> 6,
"GND"   => 3,
"GNY"   => 0,
"KRT"   => 3,
"MDT"	=> 5,

"NBB"   => 0,
"NNK"   => 0,
"NNN"   => 0,
"NNT"   => 0,
"NWY"   => 0,

"NYK"   => 1,
"RMA"   => 9,
"RNS"   => 1,

"RRS"   => 9,
"RRT"	=> 3,
"RWG"	=> 3,
"SMN"   =>  0,
"SMR"   => 3,
"SMW"	=> 5,
"SNN"   =>  0,
"SVN"   =>  0,
"SWG"   =>  2,
"TKK"	=>	3,

"TSS"	=>	3,
"VAW"   => 2,
"VHG"   => 6,
"VHK"   => 0,
"VMW"   => 4,
"VNY"	=>	0,
"VNS"   => 0,
"VYG"	=> 3,
"WHT"   => 0,
"WNY"   => 2,
"WVY"	=> 0,

# 정해진 배열, 특이적인 재결합이 잘 이루어지게 하기위해 삽입.
"GAA" => 5, # Glu
"GAT" => 1, # Asp
#"GGT" => 1, # Gly
"ATT" => 1, # Ile
#"CTG" => 1, # Leu
"TGT" => 1, # Cys
#"TGC" => 1, # Cys
"CAG" => 1, # Gln
#"AAA" => 1, # Lys
#"GAC" => 2  # Asp
);
my ($numOfMulticodons) = 0;
# set numOfMulticodons
for(keys %multicodons) {
	my ($number) = $multicodons{$_};
	$numOfMulticodons += $number;
}

my ($numOfTrial) = 5;	# repeat of random PCRed Sequence generation.
my $numOfTranslation = 5;
my (%enzymeDataBase); #DB read from file, RestructionEnzymeList, ex) "KpnI" => "GGTACC"
my $numOfEnzymeDB;

my @enzymes = ("BamHI", "NcoI", "EcoRI");	# investigate these enzyme cut sites.

# 3차 배열, 각 요소는 지정된 배열의 지정된 효소에 의한 인식부위의 첫염기
my @cutPositions = ();		# 3-dimension	[Index_Seq][Index_Enz][Index_Cut]
my @numOfCutsAllSeqs = ();	# 1-dimension	[Index_Seq]
my @pOfNoCutAllSeqs = ();	# 1-dimension	[Index_Seq]

# multibases
# elements which represent the probablilty of appearance
my @probabilityOfEachBasesOnMultibases = (
	#T, C, A, G
	[1, 0, 0, 0],	# T		index:0
	[0, 1, 0, 0],	# C		index:1
	[0, 0, 1, 0],	# A		index:2
	[0, 0, 0, 1],	# G		index:3

	[0.5, 0, 0, 0.5],	# K		index:4
	[0, 0.5, 0.5, 0],	# M
	[0.5, 0, 0.5, 0],	# W
	[0, 0.5, 0, 0.5],	# S
	[0.5, 0.5, 0, 0],	# Y
	[0, 0, 0.5, 0.5],	# R
	[0.3333, 0.3333, 0, 0.3333],	# B index:10
	[0.3333, 0, 0.3333, 0.3333],	# D
	[0.3333, 0.3333, 0.3333, 0],	# H
	[0, 0.3333, 0.3333, 0.3333],	# V
	[0.25, 0.25, 0.25, 0.25]	# N		index: 14
);

# Hash table의 요소가 배열을 지정하게 할 수 없어서 (내가 잘 몰라서),
# 일단 스트링으로 했다.
# 사용할때에는 필요에 따라서 배열로 변형 할 것 ex) @array = split(//, $basesFromMultibases{$someKey});
my %basesFromMultibases = (
   	"A" => 'A', "C" => 'C',
	"G" => 'G', "T" => 'T',
	"K" => 'GT', "M" => 'AC',
	"W" => 'AT', "S" => 'CG',
	"Y" => 'CT', "R" => 'AG',
	"B" => 'CGT', "D" => 'AGT', 
	"H" => 'ACT', "V" => 'ACG',
	"N" => 'ACGT'
);

my (%multibaseIndex) = (
	"T" => 0, "C" => 1,
	"A" => 2, "G" => 3,
	"K" => 4, "M" => 5,
	"W" => 6, "S" => 7,
	"Y" => 8, "R" => 9,
	"B" => 10, "D" => 11, 
	"H" => 12, "V" => 13,
	"N" => 14
);
my (%multibaseFromIndex) = ();
foreach (keys %multibaseIndex) {
	my $aKey = $_;
	$multibaseFromIndex{$multibaseIndex{$aKey}} = $aKey;
}

my @frequencyOfCodonInEachOfMulticodons =();

	
######## DATA
# Primers which are to be generated.
my ($PrimerFwSeq) = "";
my ($PrimerRvSeq) = "";

# PCRed sequences which are generated by random distribution of multicodons
my (@PCRedSeqs);
my ($PCRedSeq) = "";

my $input_file = '';
my $output_file = 'result.txt'; #default output file

######## Other parameters
my ($sizeOfPCRedSeq) = $numOfMulticodons * 3 + $sizeOfOverLap;
my $sizeOfFullSeqWithBothTerminals = $sizeOfPCRedSeq + (length $cloningSite5) + (length $cloningSite3);
my $sizeAdjustment = int($sizeOfFullSeqWithBothTerminals /10 /3) *3;	# 프라이머의 길이 비율을 3：2로

my ($sizeOf5End) = int($numOfMulticodons / 2) * 3;
my ($sizeOf3End) = int($numOfMulticodons / 2) * 3 + ($numOfMulticodons % 2) * 3;

# 프라이머 길이의 수정
#$sizeOf5End += $sizeAdjustment;
#$sizeOf3End -= $sizeAdjustment;

my ($sizeOfFw) = $sizeOfOverLap + $sizeOf5End;			# Number of Random Forward primer
my ($sizeOfRv) = $sizeOfOverLap + $sizeOf3End;			# Reverse primer

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
"translation"	=> "false",
"analysis" => "false",  ## 멀티코돈에서 출현 가는한 코돈들의 확률을 구한다. annealing region의 설계에 쓸 생각이다. 멀티코돈에서 자주 출현하는 코돈은 annealing region에서 사용하지 않는 것이 specific한 anneal에 중요할 것이기때문이다.
"justTestGivenFullSeq"	=> "false",
"primerGeneration"	=> "false"
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
    
    # -p : primer make from specified seq in this program
	if($temp =~ /^-\S*p\S*/){
		$flags{"primerGeneration"} = "true";
		$numOfTrial = 1;
	}
	
    # -p : primer make from specified seq in this program
	if($temp =~ /^-\S*j\S*/){
		$flags{"justTestGivenFullSeq"} = "true";
		$numOfTrial = 1;
	}
	
    # -v : visual option
	if($temp =~ /^-\S*v\S*/){
		$flags{"visual"} = "true";
	}
	# -t[Number] : show possible translations
	if($temp =~ /^-\S*t(\d*)\S*/){
		$flags{"translation"} = "true";
		$numOfTranslation = $1;
		unless (defined $numOfTranslation) { $numOfTranslation = 5; }
	}
	# -s[Number] : get specified number of sequences
	if($temp =~ /^-\S*s(\d*)\S*/){
		$numOfTrial = $1;
		unless (defined $numOfTrial) { $numOfTrial = 10; }
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
	# -a : analysis for anneal region generation
	if($temp =~ /^-\S*a\S*/){
		$flags{"analysis"} = "true";
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
	if($temp =~ /^-[vtdoajpe]*([^tdvaojpe\d])+[tdvjpaoe]*$/){
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



####################### Codon usage analysis on multi-codon, for anneal region generation
## 멀티코돈에서 출현 가는한 코돈들의 확률을 구한다. annealing region의 설계에 쓸 생각이다.
#멀티코돈에서 자주 출현하는 코돈은 annealing region에서 사용하지 않는 것이,
#specific한 anneal에 중요할 것이기때문이다.

PrintOrDisp ("Processing .... ");

my @codonUsageTable = ();    # 3-dimension array [1st base][2nd base][3rd base], uses index of %multibaseIndex

#initialize @frequencyOfCodonInEachOfMulticodons
my @base = ('T', 'C', 'A', 'G');
foreach(@base) {
	my $i = $multibaseIndex{$_};
    foreach (@base) {
    	my $j= $multibaseIndex{$_};
        foreach (@base) {
        	my $k= $multibaseIndex{$_};
            $frequencyOfCodonInEachOfMulticodons[$i][$j][$k] = 0.0;
        }
    }
}

foreach (keys %multicodons) {
    my $aMulticodon = $_;
    my $usage = $multicodons{$aMulticodon};
    if( $usage == 0 ) {
        next;
    }
    
    # 멀티코돈이 지정가능한 코돈의 종류를 나열
    my @codons = getCodonsFromMulticodon( $aMulticodon );
    
    # for debug
    if ($flags{'analysis'} eq 'true') {
        print "$aMulticodon (";
        print scalar @codons;
        print ") => ";
        &read_array(\@codons);
        print "\n";
    }
    
    # 각 코돈의 출현 확률을 @frequencyOfCodonInEachOfMulticodons 저장한다.
    my $numberOfCodons = scalar @codons;
    foreach( @codons ) {
        my $aCodon = $_;
        my @bases = split(//, $aCodon);
        
        my $index1 = $multibaseIndex{$bases[0]};
        my $index2 = $multibaseIndex{$bases[1]};
        my $index3 = $multibaseIndex{$bases[2]};
        
        $frequencyOfCodonInEachOfMulticodons[$index1][$index2][$index3] += $usage/$numberOfCodons;
    }
}

my $dummy_howMany = scalar keys %multicodons;

PrintOrDisp ("$dummy_howMany multi-codons were analyzed.\n");

# display codon usage from given multicodon usages
if( $flags{"analysis"} eq "true") {

	my $totalUsage = 0.0;
	my %sortedCodonUsage = ();
	
	PrintOrDisp("\n[CODON(Amino Acid)] : Expected Occurence\n");
	PrintOrDisp("----------------------------------------------------------------------\n");
	
	# first base
    for(my $i=0; $i<4; $i++) {
    	
    	# 3rd base
        for(my $k=0; $k<4; $k++) {
        
        	# 2nd base
            for(my $j=0; $j<4; $j++) {
            	my $aCodon = $multibaseFromIndex{$i}.$multibaseFromIndex{$j}.$multibaseFromIndex{$k};
            	my $aAA = dna2peptide($aCodon);
            	my $aCodonAAA = "[".$aCodon."(".$aAA.")]";
            	my $usage = $frequencyOfCodonInEachOfMulticodons[$i][$j][$k];
                PrintOrDisp ("$aCodonAAA %8.6f  ", $usage);
                
                # this is for sort
                if( defined $sortedCodonUsage{$usage} ) {
                	$sortedCodonUsage{$usage} .= "\t".$aCodonAAA;
                }else{
                	$sortedCodonUsage{$usage} = $aCodonAAA;
                }
                
                $totalUsage += $usage;
            }
            PrintOrDisp ("\n");
        }
        PrintOrDisp ("\n");
    }
	
	PrintOrDisp("Occurence\tCodons\n---------\t------------------------------------------------\n");
    foreach( sort numSort keys %sortedCodonUsage ) {
    	my $usage = $_;
    	my $codons = $sortedCodonUsage{$usage};
    	$codons =~ s/[\[\]]//g;
    	PrintOrDisp ("%8.6f\t$codons\n", $usage);
    }
    PrintOrDisp ("---------\t------------------------------------------------\n%8.6f\tTotal\n\n", $totalUsage);
    
    exit(1);
}


####################### Random PCRed Sequence generation
# 총코돈갯수의 10배 크기의 배열을 만들고 그 안에 랜덤으로 집어 넣는다.
# 크기를 10배로 하는 이유는 랜덤으로 추출했을 때 겹치는 offset이 생기는 확률을 낮추기 위함이며,
# 겹치는 offset을 얻으면 다시 시도하도록한다.
PrintOrDisp ("Processing .... ");
my $sizeOfDummyCodonContainer = $numOfMulticodons * 10;

srand(time);
for (my $ctr1 = 0; $ctr1 < $numOfTrial; $ctr1++) {
	
	my @dummyCodonContainer;
	my @shuffledCodons;
	my $shuffledCodonSeq = '';
	
	#fisher_yates_shuffle( \@array );    # permutes @array in place
	
	# 멀티코돈의 갯수 만큼 10배큰 배열안에 랜덤하게 분산시킨다.
	foreach (keys %multicodons) {
		for(my $usage = $multicodons{$_}; $usage >0; $usage--) {
			#push (@dummyCodonContainer, $_);
			
			my $offset;
			do {
				$offset = int(rand($sizeOfDummyCodonContainer));	#generate random offset
			} while (defined $dummyCodonContainer[$offset]);		# get another offset if not empty element
			
			$dummyCodonContainer[$offset] = $_;
		}
	}

	# 정의 되지 않은 빈 요소를 제외한 랜덤하게 분산된 요소들만 모아서 새로운 배열을 만든다.
	for(my $ctr2 = 0; $ctr2 < $sizeOfDummyCodonContainer; $ctr2++) {
		my ($aCodon) = $dummyCodonContainer[$ctr2];
		if (defined $aCodon) {
			push (@shuffledCodons, $aCodon);
		}
	}
	
	# 연속된 코돈이 오지 않는 것이 바람직한가?
	# 나중에 제한효소 부위가 적게나오는 셀렉션을 할 계획인데 그러면, 어떠한 반복패턴이 유리하게 추출될 가능성이 있다.
	# 이중에서 연속된 코돈 특히 경우의수가 작은 멀티코돈이 연속되게 나올 가능성이 높을 것 같다.
	# 이걸 줄이기 위해 다시 한번 섞는다.
	# 단, 출현빈도가 적은 코돈만 다시 뒤섞는다.
	my $lowFrequencyCodon = $numOfMulticodons * 0.04;
	for(my $ctr2=0; $ctr2<@shuffledCodons - 1; $ctr2++) {
		my $currentCodon = $shuffledCodons[$ctr2];
		
		# 코돈의 출현 빈도가 낮은 코돈이거나, 멀티코돈이 아닌 코돈(A,C,G,T로만 지정된 코돈)일때, +세번째 염기만 두가지 염기가 될 가능성을 가질때
		if ($multicodons{$currentCodon} < $lowFrequencyCodon || $currentCodon =~ /[ACTG][ACTG][ACTGSWRYKM]/ ) {
			#print "$currentCodon\n";
			
			# check if this low codon followed by same codon
			if($currentCodon eq $shuffledCodons[$ctr2 + 1]) {
				#print "$shuffledCodons[$ctr2] shuffled\n";
				#exchange with randomly chosen codon
				my $indexRandom = int(rand(@shuffledCodons));
				($shuffledCodons[$ctr2], $shuffledCodons[$indexRandom]) = ($shuffledCodons[$indexRandom], $shuffledCodons[$ctr2]);
			}
		}else {
			#print " High frequency : $currentCodon\n";
		}
	}
	
	# 랜덤소팅된 배열의 요소를 뭉쳐서 DNA를 작성한다.
	$shuffledCodonSeq = join '', @shuffledCodons;		# generate random seq from multi-codons
	
	# 중간에 링커시퀀스를 끼워넣는다.
	$PCRedSeq = substr ($shuffledCodonSeq, 0, $sizeOf5End);	# set 5' seq  
	$PCRedSeq .= $LinkerSeq;						# concatanate Linker with 5' seq
	$PCRedSeq .= substr($shuffledCodonSeq, $sizeOf5End);		# concatanate 3' seq
	
	push (@PCRedSeqs, $PCRedSeq);
}
$dummy_howMany = scalar @PCRedSeqs;
PrintOrDisp ("$dummy_howMany Sequences were generated.\n");

if ($flags{'justTestGivenFullSeq'} eq 'true') {
	$SequenceForjustTestGivenFullSeq =~ s/\s//g;
	$PCRedSeqs[0] = $SequenceForjustTestGivenFullSeq;
}

####################### Primer generation
# 작성한 랜덤 시퀀스를 중합할 수 있은 프라이머를 생성한다.
# 이때 지정된 Annealing region이 랜덤 시퀀스의 중간에 오게하며,
# 각 프라이머의 3'에 Annealing region이 위치,
# 5'에 지정된 Restruction site를 붙인다.
PrintOrDisp ("Processing .... ");

# do it for every seqs
@primersOfFw = ();
@primersOfRv = ();
for(my $i=0; $i< @PCRedSeqs; $i++ ) {
	$PCRedSeq = $PCRedSeqs[$i];
	
	my $primerFw;
	my $primerRv;
	my $sizeOfFw;
	my $sizeOfRv;
	
	if( $PCRedSeq =~ m/$LinkerSeq/ ) {
	
		# Add Annealing Region on Random Region
	 	$primerFw = $` . $&;	# 5' + Annealing region
		$primerRv = $& . $';	# Annealing region + 3'
		
		# Add Cloning Site
		$primerFw = $cloningSite5 . $primerFw;
		$primerRv = $primerRv . $cloningSite3;
		
		# get reverse-complement seq
		$primerRv = GetReverseComplement($primerRv);
	
		# Save on primer database
		$primersOfFw[$i] = $primerFw;
		$primersOfRv[$i] = $primerRv;
		
		$sizeOfFw = length $primerFw;
		$sizeOfRv = length $primerRv;
	}else{
		print "		Fatal Error	: Cannot find Annealing seq. [Ox01]\n";
		exit();
	}
}
$dummy_howMany = scalar @PCRedSeqs;
PrintOrDisp ("$dummy_howMany primer pairs were generated.\n");


####################### Possible Annealing Region Check
# PCR을 이용하여 프라이머로 Random Sequence를 중합할때 생길 수 있는 부산물의 종류와
# 확률을 계산한다.
# Annealing Region (Linker Seq)가 Random Region에 출현하면 Annealing 시 부적절한 중합이 일어나고,
# 반복적인 PCR을 할 경우(cycle 횟수가 2이상일경우), 첫번째 산물들끼리 중합하여 더 큰 부산물이 생성될 수 있다.
# Restruction Site 체크와 같은 알고리즘을 사용한다.
# 배열체크 알고리즘은 서브루팅으로 돌리자.

PrintOrDisp ("Processing .... ");

# set anneal seqs to search
%annealSeqs = ();
# Reverse primer
$annealSeqs{'reverse'} = $LinkerSeq;
# Forward primer
$annealSeqs{'forward'} = reverse $LinkerSeq;

#my $direction = 'forward'; # or 'reverse'

my $annealMaxLength = length $LinkerSeq;    # 8 codons
my $annealMinLength = 7; #int ($annealMaxLength / 2);


# check for all PCRedSeqs
@annealPositions = ();
for(my $i=0; $i< @PCRedSeqs; $i++ ) {
	
	# set searching region on
	my $fullPCRedSeq = $cloningSite5 . $PCRedSeq . $cloningSite3;
	my $sizeOfFullPCRedSeq = length $fullPCRedSeq;
	
	if ( $flags{'debug'} eq 'true' ) {
		my $seqByThreeLetters = GetSeqByThreeLetters($fullPCRedSeq);
		PrintOrDisp ("$seqByThreeLetters\n");
	}
	
	# check from here
	my $numOfAnneals = 0;
	$pOfNoMatchAllSeqs[$i] = 1.0;	# init for calc.
	$numOfMatchesAllSeqs[$i] = 0;	# init for calc.
	foreach (keys %annealSeqs) {
		my $direction = $_;
		my $j = 0;
		if ($direction eq 'forward') {
		    $j = 1;
		}else{
		    $j = 0; # reverse
		}
		
		# start from minimal matching
		#my $recognitionSite = substr ($annealSeqs[$j], 0, $annealMinLength);
		
		my $probability = 1.0;			# p of Cut
		my $pOfNoMatchAtSinglePosition = 1.0;
		my $pOfNoMatches = 1.0;			# init for calc.
		my $numOfMatches = 0;

		# KMP 방법을 쓰면 더 빠를 수도 있겠지만, 4개밖에 없는 DNA배열특성상 효율이 크게 향상되지는 않을 것 같다.
		# 그냥 알기 쉽게 차근차근 비교하는 방법으로 search하자.
		
		my @recognitionSiteInArray = split (//, $annealSeqs{$direction});
		my @PCRedSeqInArray = '';
		
		if ($direction eq 'reverse') {
		    @PCRedSeqInArray = split (//, $fullPCRedSeq);
		}elsif ($direction eq 'forward') {
		    # 역방향으로 비교해야하므로 비교할 대상을 뒤집어 놓는다.
		    @PCRedSeqInArray = split (//, reverse $fullPCRedSeq);
		    
		}

		# 배열의 처음부터 끝까지 비교
		for (my $ctr1 = 0; $ctr1 < $sizeOfFullPCRedSeq - $annealMinLength; $ctr1++) {
		
		    # for comparison
		    my $DNA1 = '';
		    my $DNA2 = '';
		    my $p = 0;
		    
			# 지정된 효소의 인식부위와 일대일비교
			$probability = 1.0;
			for (my $ctr2 = 0; $ctr2 < $annealMaxLength; $ctr2++) {
				
				# 참조범위를 넘어가면 루프를 끝낸다. @PCRedSeqInArray의 마지막부분때문에 이렇게 처리한다.
				# 전체크기 - $annealMinLength 까지 하지만, 최대 $annealMaxLength까지 비교하기에,
				# PCRedSeqInArray의 크기를 초과해서 비교할 수 도 있기때문이다.
				if($ctr1 + $ctr2 >= $sizeOfFullPCRedSeq) {
				    $ctr2 = $annealMaxLength;   # for ending the loop
				    $probability = 0;
				}else{
				    $DNA1 = $PCRedSeqInArray[$ctr1 + $ctr2];
				    $DNA2 = $recognitionSiteInArray[$ctr2];
				
				    if(defined $DNA1 && defined $DNA2) {
				        $p = GetPOfSameBaseOccurance($DNA1, $DNA2); # single base
				    }else{
				        if(defined $DNA1) {
				            print "DNA2 is not defined. i=$i direction=$direction ctr1=$ctr1 ctr2=$ctr2\n";
				        }else{
				            print "DNA1 is not defined. i=$i direction=$direction ctr1=$ctr1 ctr2=$ctr2\n";
				        }
				    }
				    $probability *= $p;	# get probility of match
				}
				
				if($probability == 0) {
					# No match
					$ctr2 = $annealMaxLength;	# End the roop
				
				# $probability > 0
				}else{
					if( $ctr2 == $annealMinLength-1 ) {	# only increase once for each matches
						
						# cacl. probability that enzyme cut never occurs
						$pOfNoMatchAtSinglePosition = (1.0 - $probability);
						$pOfNoMatches *= $pOfNoMatchAtSinglePosition;
							
						$annealPossibility[$i][$j][$numOfMatches] = $probability;
							
						# record cut positions at three-dim-array[Index_Seq][Index_Enz][Index_Anneal]
						$annealPositions[$i][$j][$numOfMatches] = $ctr1;
					
					    # record matched length
					    $annealLengthes[$i][$j][$numOfMatches] = $ctr2+1;
						
						# increase cut count for next
						$numOfMatches++;
							
					}elsif( $ctr2 > $annealMinLength ) {
						$annealPossibility[$i][$j][$numOfMatches-1] += $probability;
							
						# record cut positions at three-dim-array[Index_Seq][Index_Enz][Index_Anneal]
						$annealPositions[$i][$j][$numOfMatches-1] = $ctr1;
					
						# record matched length
						$annealLengthes[$i][$j][$numOfMatches-1] = $ctr2+1;
					}
				} # end of if-else()
			
			} # end of for($ctr2)
			
		} # end of for($ctr1)
		
		$pOfNoMatchAllSeqs[$i] *= $pOfNoMatches;
		$numOfMatchesAllSeqs[$i] += $numOfMatches;

	} # end of searching each primer matching
	
	# for debugging
	if ( $flags{'debug'} eq 'true') {
		&read_array(\@annealLengthes); # show the whole array
	}
	
} # end of searching each seq

PrintOrDisp ("Possible annealing positions are calculated.\n");

####################### Restruction Enzyme Site Check
# 작성한 랜덤소팅된 멀티코돈 배열에서 지정한 제한효소 부위가 등장하는지 검사한다.
# 부위의 갯수가 없는 배열조합을 선택출력한다.
PrintOrDisp ("Processing .... ");

# check for all PCRedSeqs
@cutPositions = ();
my %cutPatternsInAllEnz;

for(my $i=0; $i< @PCRedSeqs; $i++ ) {
	$PCRedSeq = $PCRedSeqs[$i];
	
	if ( $flags{'debug'} eq 'true' ) {
		my $seqByThreeLetters = GetSeqByThreeLetters($PCRedSeq);
		PrintOrDisp ("$seqByThreeLetters\n");
	}
	
	# check for all specified restruction enzymes
	
	my $cutNumberOfAllEnzyme = 0;
	$pOfNoCutAllSeqs[$i] = 1.0;	# init for calc.
	$numOfCutsAllSeqs[$i] = 0;	# init for calc.
	for(my $j=0; $j<@enzymes; $j++) {
		my ($enzyme) = $enzymes[$j];
		my ($recognitionSite) = $enzymeDataBase{$enzyme};
		my $sizeOfRecognitionSite = length $recognitionSite;
		
		my $probability = 1.0;			# p of Cut
		my $pOfNoCutAtSinglePosition = 1.0;
		my $pOfNoCuts = 1.0;			# init for calc.
		my $numOfCuts = 0;
		
		#%cutPatternsForEachEnz = ();
		
		if ( $flags{'visual'} eq 'true' ) {
			PrintOrDisp ("\n#%03d \n", $i+1);
			PrintOrDisp ("$enzyme: \t$recognitionSite\n");
		}
		
		# KMP 방법을 쓰면 더 빠를 수도 있겠지만, 4개밖에 없는 DNA배열특성상 효율이 크게 향상되지는 않을 것 같다.
		# 그냥 알기 쉽게 차근차근 비교하는 방법으로 search하자.
		my (@PCRedSeqInArray) = split (//, $PCRedSeq);
		my (@recognitionSiteInArray) = split (//, $recognitionSite);
		
		
		#### 여기서부터 서브루틴으로 넘겨
		# position[] = findAllMatches(String1, String2)
		#@position = findAllMatches($PCRedSeq, $recognitionSite);
####################### SUB FindAllMatches - 이런 이런거 만들 시간이 없다.
# 그냥 보기 안 좋더라도 메인 함수에서 돌리자. 나주에 시간나면 만들자 ㅠㅠ
# Sub-function
# FindAllMatches
# March 9th, 2012 Yusung KIM
# 첫째 배열안에서 둘째배열이 출현하는지 또는 출현 할 수 있는지 계산하여 그 위치를 배열로 전달한다.
# usage: @position = findAllMatches($PCRedSeq, $recognitionSite);
#sub FindAllMatches {
#}

		# 배열의 처음부터 끝까지 비교
		my $startFrom = 0;
		my $endTo = $sizeOfPCRedSeq - $sizeOfRecognitionSite;
		if($flags{'justTestGivenFullSeq'} eq 'true') {
			$endTo = (length $PCRedSeqs[0]) - $sizeOfRecognitionSite;
		}
		
		for (my $ctr1 = $startFrom; $ctr1 < $endTo; $ctr1++) {
		
			if ( $flags{'debug'} eq 'true' ) { print "[pos:$ctr1] "; } # for debugging
			
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

				# record recognition pattern
				my $aPattern;
				for( $ctr1 .. $ctr1+$sizeOfRecognitionSite-1 ) {
					$aPattern .= $PCRedSeqInArray[$_];
				}
				#$cutPatternsInAllEnz{substr($PCRedSeq, $ctr1, $sizeOfRecognitionSite)} = $probability;
				$cutPatternsInAllEnz{$aPattern} = $probability;
				
				# record cut positions at three-dim-array[Index_Seq][Index_Enz][Index_Cut]
				$cutPositions[$i][$j][$numOfCuts] = $ctr1;
				
				# increase cut count
				$numOfCuts++;
				
				# cacl. probability that enzyme cut never occurs
				$pOfNoCutAtSinglePosition = (1.0 - $probability);
				$pOfNoCuts *= $pOfNoCutAtSinglePosition;
				
				if ( $flags{'visual'} eq 'true' ) {
					my $seq = substr ( $PCRedSeq, $ctr1, $sizeOfRecognitionSite );
					PrintOrDisp( "\tcut at %03d\t", $ctr1);
					PrintOrDisp( "\"$seq\"\t");
					PrintOrDisp( "p=%.6f\t", $probability);
					PrintOrDisp( "pOfNoCut=%.4f\n", $pOfNoCuts);
				} # end of if
				
			} # end of if
			
		} # end of for()
		
		#$cutPatternsForEnzymeOf{$enzyme} = \%cutPatternsForEachEnz;
		#push (@cutPatternsForAllEnz, \%cutPatternsForEachEnz);
		
		$pOfNoCutAllSeqs[$i] *= $pOfNoCuts;
		$numOfCutsAllSeqs[$i] += $numOfCuts;

		my $noCutPercentage = $pOfNoCuts * 100;
		if ( $flags{'visual'} eq 'true' ) {
			PrintOrDisp (" => Possible cuts of %8s", $enzymes[$j]);
			PrintOrDisp (": $numOfCuts\tpNoCut: $noCutPercentage %%\n");
		}

	} # end of searching each Enzyme
	
	#@cutPositions = sort numSort @cutPositions;
	
	# for debugging
	if ( $flags{'debug'} eq 'true') {
		&read_array(\@cutPositions); # show the whole array
	}
	
} # end of searching each seq

$dummy_howMany = scalar @PCRedSeqs;
my $dummy_howMany2 = scalar @enzymes;
PrintOrDisp ("$dummy_howMany2 enzymes (");
for(my $i=0; $i<$dummy_howMany2-1; $i++) { PrintOrDisp ("$enzymes[$i], "); }
PrintOrDisp ("$enzymes[$dummy_howMany2-1]) were tested on each of $dummy_howMany sequence(s).\n");

######################## Generate all codons from a given multicodon
# sub get codons from given single multicodon
sub getCodonsFromMulticodon {
    my ($aMulticodon) = @_;
    my @codons = ();       # save all codons here
    my $aCodon = '';        # for calc.
    
    # get multibases on each position of codon
    # 만약 멀티코돈이 정의 되지 않았거나 길이가 짧으면 에러를 발생.
    if((!defined $aMulticodon) || (length $aMulticodon) != 3) {
        print "Warning args: sub getCodonsFromMulticodon\n";
        return;
    }
    my $thirdBase = substr($aMulticodon, 2, 1);
    my $secondBase = substr($aMulticodon, 1, 1);
    my $firstBase = substr($aMulticodon, 0, 1);

    my @firstBases = getBasesFromMultibase ($firstBase);
    my @secondBases = getBasesFromMultibase ($secondBase);
    my @thirdBases = getBasesFromMultibase ($thirdBase);
    
    # make codons by combination
    for(my $ctr1=0; $ctr1<@firstBases; $ctr1++) {
        for(my $ctr2=0; $ctr2<@secondBases; $ctr2++) {
            for(my $ctr3=0; $ctr3<@thirdBases; $ctr3++) {
                # generate a codon by concatenating each base
                $aCodon = $firstBases[$ctr1] . $secondBases[$ctr2] . $thirdBases[$ctr3];
                
                # save it
                push (@codons, $aCodon);
                #if ($aCodon eq 'TAA') {
                #	print "____ $aMulticodon ____\n";
                #}
            }
        }
    }
    
    #&read_array(\@codons);
    return @codons;
}

######################## Generate all bases from a given multibase
# sub get codons from given single multibase
sub getBasesFromMultibase {
    my ($aMultibase) = @_;
    
    my @allBases = split(//, $basesFromMultibases{$aMultibase});
    
    #&read_array(\@allBases);
    return @allBases;
}

######################## Sort for Selection of best sequence
# sort seqs by High pOfNoCut
# sort algorithm : insertion sort

# calc. num of Frame shift cuts
my @numOfFrameShiftCutsAllSeq = ();
for(my $i=0; $i<$numOfTrial; $i++) {
	my $numOfFrameShiftCuts = 0;
	for (my $j = 0; $j < @enzymes; $j++) {
		for (my $k = 0; $k < $#{$cutPositions[$i][$j]} + 1; $k++) {
			if ($cutPositions[$i][$j][$k] % 3 != 0) {
				$numOfFrameShiftCuts++;
			}
		}
	}
	$numOfFrameShiftCutsAllSeq[$i] = $numOfFrameShiftCuts;
}

for(my $i=0; $i<$numOfTrial; $i++) {
	my $current = 0.0;
	my $max = 0.0;
	my $indexMax = 0;
	
	
	$indexMax = $i;
	$max = $numOfFrameShiftCutsAllSeq[$i];
	
	#search min
	for(my $j=$i+1; $j<$numOfTrial; $j++) {
		
		$current = $numOfFrameShiftCutsAllSeq[$j];
		
		# set new min
		if($current > $max) {
			$max = $current;
			$indexMax = $j;
		}
	}
	
	# swap min with new min
	if( $i != $indexMax ) {
		# Swap Sequences
		($PCRedSeqs[$i], $PCRedSeqs[$indexMax]) = ($PCRedSeqs[$indexMax], $PCRedSeqs[$i]);
		
		# Swap Total Cuts
		($numOfCutsAllSeqs[$i], $numOfCutsAllSeqs[$indexMax]) = ($numOfCutsAllSeqs[$indexMax], $numOfCutsAllSeqs[$i]);
	
		# Swap Total Frame shift cuts
		($numOfFrameShiftCutsAllSeq[$i], $numOfFrameShiftCutsAllSeq[$indexMax]) = ($numOfFrameShiftCutsAllSeq[$indexMax], $numOfFrameShiftCutsAllSeq[$i]);
		
		# Swap Cut Position
		($cutPositions[$i], $cutPositions[$indexMax]) = ($cutPositions[$indexMax], $cutPositions[$i]);
	
		# Swap probability
		($pOfNoCutAllSeqs[$i], $pOfNoCutAllSeqs[$indexMax]) = ($pOfNoCutAllSeqs[$indexMax], $pOfNoCutAllSeqs[$i]);
	}
}


#################################### Report
PrintOrDisp( "\n<Report>\n" );
PrintOrDisp( "=====================================================================\n" );
my $seqByThreeLetters = '';
my @cutScreens = ();
for(my $i=0; $i<$numOfTrial; $i++) {

	PrintOrDisp ("#%03d ", $i+1);
	PrintOrDisp ("(pNoCut=%.2f%%)\t", $pOfNoCutAllSeqs[$i]*100);
	
	# print cut statistics
	my $numOfEnzymes = $#enzymes + 1;	#get size of an array, or you can do 'scalar @enzymes'
	# $numOfEnzymes = @enzmyes;
	# $numOfEnzymes = scalar @enzymes;
	
	for(my $j=0; $j<@enzymes; $j++) {
		my $enzymeNum = $j+1;
		
		my $cuts = $#{$cutPositions[$i][$j]} + 1;
		#my $cuts = $#{$cutPositions[$i]} + 1;
		
		# ex) #1 NdeI (6), #2 EcoRI (5), #3 ...
		PrintOrDisp ("#$enzymeNum:$enzymes[$j]($cuts)");
		
		if($j<$numOfEnzymes-1) {
			PrintOrDisp (" + ");
		}
	}
	
	# print total cuts
	PrintOrDisp ("\t Total %d Cuts", $numOfCutsAllSeqs[$i]);
	
	# print total cuts which induce frame shift
	PrintOrDisp ("(%d Frame Shift Cuts)\n", $numOfFrameShiftCutsAllSeq[$i]);
	
	# print sequences
	# -j 옵션에도 대응을 해주기 위해 이렇게 변경함.
	if($flags{'justTestGivenFullSeq'} eq 'true') {
		$sizeOfPCRedSeq = length $PCRedSeqs[0];
	}
	
	PrintOrDisp ("          [DNA Sequence] $sizeOfPCRedSeq base\n");
	$seqByThreeLetters = GetSeqByThreeLetters($PCRedSeqs[$i]);
	PrintOrDisp ("          $seqByThreeLetters\n");
	
	for (my $j = 0; $j < @enzymes; $j++) {
		my $cutScreen = '';		# cut position will be set on this string
		my $nextPosition;
		my $indexNext = 0;		# index for @totalCutPositions
	
		# set first cutPosition
		if( defined $cutPositions[$i][$j][0]) {
			$nextPosition = $cutPositions[$i][$j][0];
		}else{
			$nextPosition = $sizeOfPCRedSeq;	# if no cut position
		}
		
		# print every cut position on the top of seq
		for (my $k = 0; $k< $sizeOfPCRedSeq; $k++) {
			if($k<$nextPosition) {
				$cutScreen .= "-";
			}elsif($k == $nextPosition){
				my $EnzymeSite = $enzymeDataBase{$enzymes[$j]};
				
				if(++$indexNext < $#{$cutPositions[$i][$j]}+1) {
					$nextPosition = $cutPositions[$i][$j][$indexNext];
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
		PrintOrDisp ("$seqByThreeLetters\n");
	}	
		
	################## print possible anneals
	# print sequences
	my $fullPCRedSeq = $cloningSite5 . $PCRedSeqs[$i] . $cloningSite3;
	my $sizeOfFullPCRedSeq = length $fullPCRedSeq;
	my $frameShift = (length $cloningSite5) % 3;
	
	PrintOrDisp ("[PCRed Seq] $sizeOfFullPCRedSeq base  \tDesignated Anneal Seq.: 5'-$LinkerSeq-3'\n");
	$seqByThreeLetters = GetSeqByThreeLetters($fullPCRedSeq, $frameShift);
	PrintOrDisp ("          $seqByThreeLetters\n");
	
	foreach (keys %annealSeqs) {
	    my $direction = $_;
		my $matchScreen = '';		# cut position will be set on this string
		my $nextPosition;
		my $indexNext = 0;		# index for @totalCutPositions
		
		my $j = 0; # for anneal direction
		if($direction eq 'forward') {
		    $j = 1;
		}else{
		    $j = 0; # reverse
		}
	
		my $annealDirectionForPrint = "[" . $direction . "]";
		PrintOrDisp ("%9s ", $annealDirectionForPrint);
		
		# set first cutPosition
		if( defined $annealPositions[$i][$j][0]) {
			$nextPosition = $annealPositions[$i][$j][0];
		}else{
			PrintOrDisp ("          NO Possible Anneal !\n");
			next;
		}
		
		# print every cut position on the top of seq
		for (my $k = 0; $k< $sizeOfFullPCRedSeq; $k++) {
			if($k<$nextPosition) {
				$matchScreen .= "-";
			}elsif($k == $nextPosition){
			    my $sizeOfMatchedRegion = $annealLengthes[$i][$j][$indexNext];
				my $matchedRegion = substr($annealSeqs{$direction}, 0, $sizeOfMatchedRegion);
				#print "pos:$k". "\tsize:$sizeOfMatchedRegion". "\tSeq:$matchedRegion". " \n";
				
				# if there is more anneals
				if(++$indexNext < $#{$annealPositions[$i][$j]}+1) { # (current index + 1) < (size of array)
					$nextPosition = $annealPositions[$i][$j][$indexNext];
					if($k <= $nextPosition - $sizeOfMatchedRegion ) {
						$matchScreen .= $matchedRegion;
						$k = $k - 1 + $sizeOfMatchedRegion;
					}else{
						$matchScreen .= "|";	# if another site exist closly behind
					}
					
				# if this is the last anneal
				}else{
					$matchScreen .= $matchedRegion;
					$nextPosition = $sizeOfFullPCRedSeq;	# no more cut, end the loop
					$k = $k - 1 + $sizeOfMatchedRegion;
				}
			}
		}
		
		# print cut position by three letter spaced way
		if($direction eq 'forward') {
		    $matchScreen = reverse $matchScreen;    # reverse the seq because I calc.ed this reversly.
		}
		my $frameShift = (length $cloningSite5) % 3;
		$seqByThreeLetters = GetSeqByThreeLetters($matchScreen, $frameShift);
		PrintOrDisp ("$seqByThreeLetters\n");
	}
	PrintOrDisp("\n");


	# make primers
	if ( $flags{'primerGeneration'} eq 'true') {
		my $primerFw;
		my $primerRv;
		my $sizeOfFw;
		my $sizeOfRv;
	
		my $fullSequenceWithBothTerminals = $cloningSite5. $PCRedSeqs[$i] . $cloningSite3;
		if( $fullSequenceWithBothTerminals =~ m/$LinkerSeq/ ) {
	 		$primerFw = $` . $&;	# 5' + Annealing region
			$primerRv = $& . $';	# Annealing region + 3'
			$sizeOfFw = length $primerFw;
			$sizeOfRv = length $primerRv;
		}else{
			print "		Fatal Error	: Cannot find Annealing seq. [Ox02]\n";
		}
	
		# get reverse-complement seq
		$primerRv = GetReverseComplement($primerRv);
	
		PrintOrDisp ( "[Primer Fw] %d bases\n", $sizeOfFw );
		PrintOrDisp ( "%s\n", GetSeqByThreeLetters($primerFw) );
		PrintOrDisp ( "[Primer Rv] %d bases\n", $sizeOfRv);
		PrintOrDisp ( "%s\n\n", GetSeqByThreeLetters($primerRv) );
}

	# debugging
	#if($flags{'debug'} eq 'true') {
	#	for (my $j=0; $j<@totalCutPositions; $j++) {
	#		print " $totalCutPositions[$j]";
	#	}
	#	print "\n";
	#}
	
	# display all sequences if visual option is on
	if($flags{'translation'} eq 'true') {
		PrintOrDisp ("\nPossible Translations...\n");
		if( defined $fullPCRedSeq ) {
			TestTranslation ($fullPCRedSeq, $numOfTranslation);
		}else{
			TestTranslation ($PCRedSeqs[$i], $numOfTranslation);
		}
	}
	PrintOrDisp ("\n");
	
}

# make primers for justTestGivenFullSeq
if ( $flags{'justTestGivenFullSeq'} eq 'true') {
	my $primerFw;
	my $primerRv;
	my $sizeOfFw;
	my $sizeOfRv;
	
	my $fullSequenceWithBothTerminals = $cloningSite5. $PCRedSeqs[0] . $cloningSite3;
	
	if( $fullSequenceWithBothTerminals =~ m/$LinkerSeq/ ) {
	 	$primerFw = $` . $&;	# 5' + Annealing region
		$primerRv = $& . $';	# Annealing region + 3'
		$sizeOfFw = length $primerFw;
		$sizeOfRv = length $primerRv;
	}else{
		print "		Fatal Error	: Cannot find Annealing seq. [Ox02]\n";
	}
	
	# get reverse-complement seq
	$primerRv = GetReverseComplement($primerRv);
	
	PrintOrDisp ( "[Primer Fw] %d bases\n", $sizeOfFw );
	PrintOrDisp ( "%s\n", GetSeqByThreeLetters($primerFw) );
	PrintOrDisp ( "[Primer Rv] %d bases\n", $sizeOfRv);
	PrintOrDisp ( "%s\n", GetSeqByThreeLetters($primerRv) );
}

# Report possible Enzyme cut pattern 
PrintOrDisp("\nPossible cut seqs and occurence for ");
for (my $i=0; $i < @enzymes; $i++) {
	my $anEnzyme = $enzymes[$i];
	
	#my $patternsForEachEnz = $cutPatternsForAllEnz[$i];
	#my $patternsForEachEnz = \%{$cutPatternsForEnzymeOf{$anEnzyme}};
	
	PrintOrDisp( "[$anEnzyme]: $enzymeDataBase{$anEnzyme}   ");
	
	#foreach (keys %{$patternsForEachEnz}) {
	#	PrintOrDisp( "\t$_ : ${$patternsForEachEnz}{$_}\n");
	#}
}
PrintOrDisp("\n");
foreach (keys %cutPatternsInAllEnz) {
	PrintOrDisp( "%9s : $cutPatternsInAllEnz{$_}\n", $_);
}
PrintOrDisp("\n");
		
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
	my $frameShift = shift;	# frame shifted, ex) x number of bases are added 
	my $readingFrame = 0;   # reading frame
	my $length = length $seq;
	
	if ($length < 1) {
		return;
	}
	
	# frame adjustment by specified reading frame
	if (defined $frameShift) {
	    #calc. frame to the right frame
		$readingFrame = (3-$frameShift) % 3;
		
		if ($readingFrame == 1) {
			$seq = ' '.$seq;    #add one space
		}elsif( $readingFrame ==2 ) {
			$seq = '  '.$seq;   #add two spaces
		}
		
	}else {
		$readingFrame = 0;
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

	my $probability = $probabilityOfEachBasesOnMultibases[$indexMultibase][$indexBase];
	
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

# usage : &read_array(\@annealLengthes); # show the whole array
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
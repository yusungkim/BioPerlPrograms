####### RandomDNA2Peptide.pl
####### 9.Apr.2013
####### usage: RandomDNA2Peptide [options] [-o Output_File] [-e Enzyme]
####### 하는일 : 입력된 Multi-codon sequence를 이용해 지정된 갯수만큼의 DNA를 생성한 후, Peptide를 출력한다.


sub version {
    print "\n";
    print "#############################################################\n";
    print "##         RandomDNA2Peptide  Version 1.0                  ##\n";
    print "##     .............................................       ##\n";
    print "##             All Copywrights are reserved by             ##\n";
    print "##           Yusung KIM, 2013.04.09 ~ 2013.04.10           ##\n";
    print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com    ##\n";
    print "#############################################################\n";
}
sub usage{
    print "\nusage: perl $0 [Input_File] [-o Output_File] [options]\n";
    print "ex)\n perl $0 INPUT.txt -o RESULT_Translated.txt -t5\n";
    print "options:\n";
	print "	-t	-> specify trial number of possible translation\n";
	print "        ex) -t5\n";
    print "	-o	-> specify output file to save the results\n";
    print "		   ex) -o output.txt\n";
    print "	-c 	-> display generated dna sequences used for translations\n";
    print "	-V	-> display version\n";
    print "	-h/?	-> display this\n";
}


################################### MAIN #################################

#!/usr/bin/perl -w
use warnings;
use strict;

######## DATA Files
my $input_file = 'INPUT_Multi-codon_Sequence.txt'; #default input file
my $output_file = 'RESULT_Translated.txt'; #default output file


# parameters
my $numOfNucleotides = 0;
my $numOfCodons	= 0;
my $numOfTranslation = 5;	# default trials of translation


# default input multi-dna sequence

# 20130723 pUC-Nark12(14)+RND
my $MultiCodonSequence = "ATGGACGGGAATTCAGACVHKTATRNSCAGTGGTTGVNHVNHRRSGGACCGVNHTCTGGTCGTCCACCACCARNSGGARRSKCTRNSGGTRRTVHKTWTRNSVAWTGGMWGVNHVNHRRSGGTCCAVNHTCTGGACGTCCTCCTCCGTCTGGTTCCGGARNSGGARRTVHKTATRNSCAATGGCTCVNHVNHRRSGGACCGVNHAGTGGTVRWCCTCCGCCARNSGGTRNSKCTRRSGGTRRTVHKTWYRNSVAWTGGTTAVNHVNHRRSGGACCGVNHAGCGGACGTCCACCTCCTRNSGGTRRSGGTGCTGCTGGATCCGGTGGTTCATCTGGAGCTAGCAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCTACAATCGGAAAACTCACCCTTAAATTTATTTGCACTACTGGAAAACTGCCTGTTCCGTGGCCAACACTTGTCACTACTCTGACCTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGTCATGACTTTTTCAAGAGTGCCATGCCGGAAGGTTATGTACAGGAACGCACTATTTCTTTCAAAGATGACGGGAAATACAAGACGCGTGCTGTAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAGGGTACTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAATACAACTTTAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCACAGTTCGCCACAACGTTGAAGATGGTTCCGTTCAACTTGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCTACACAAACTGTCCTTTCGAAAGATCCGAACGAAAAGCGGGACCACATGGTGCTGCACGAGTACGTGAACGCCGCTGGAATCACACATGGTGGTACCGGTGAAAACCTGTATTTTCAGGGCTGCAATGACATGACTCCAGAGCAAATGGCTACAAATGTGAACTGTTCCAGCCCTGAGCGTCACACACGCAGTTATGATTACATGGAAGGTGGGGATATTCGTGTGCGTCGTCTCTTCTGTCGCACACAGTGGTATCTGCGTATTGATAAACGTGGCAAAGTAAAAGGGACCCAAGAGATGAAGAATAATTACAATATCATGGAAATCCGTACAGTGGCAGTTGGAATTGTGGCAATCAAAGGGGTGGAAAGTGAGTTCTATCTTGCAATGAACAAGGAAGGAAAACTCTATGCAAAGAAAGAATGCAATGAAGATTGTAACTTCAAAGAACTGATTCTGGAAAACCATTACAACACCTATGCATCAGCTAAATGGACACACAACGGAGGGGAAATGTTTGTTGCCTTAAATCAAAAGGGGATTCCTGTACGTGGAAAAAAAACGAAGAAAGAACAAAAAACAGCCCACTTTCTTCCTATGGCAATCACTTAA";

#my $MultiCodonSequence ="ATCCAAGAAGGAGATSMWVYGRRSWNYRRSGAWRWGVYGVHGSMRRWGMDTRRTVMWGAADHTVHGVHGMDTSMWRRTVHGVMWTGTDHTRMARRTVMWRMASMWSMWTKKRMASWGRRSRMANYKRRSDHTKRTTSSGNDTKKSWGRMAGAWTSSCTGTCTAAAGACCCTAATGAGAAACGCGACRWGWNYKRTMDTGAWDHTDHTRMAMDTVYGGAARRSVMWVAWGNDGAAGATRNSGNDGAWRRSRMARRSSMRRRSRMAGAWDHTTGGGAWTKKVAWCAGRMAMDTVHGTSSGAARRSDHTSMRGAAATTDHTKRTGAACTTGTATCTGTCGAATTC";
$MultiCodonSequence =~ s/\s//g;					# Remove spaces
$MultiCodonSequence =~ s/[^A-Za-z]//g;			# Remove non-alphabet
$MultiCodonSequence =~ tr/[a-z]/[A-Z]/;			# To Upper Character
#$MultiCodonSequence =~ s/[^AGCTNYRSWHVDB]//ig;	# Remove everything except DNA

$numOfNucleotides = length $MultiCodonSequence;
$numOfCodons = scalar ($numOfNucleotides / 3);

##################### Parameters privatly created by me
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


	
####### CHECKING OPTIONS

#if( @ARGV == 0 ) {
#	version();
#	usage();
#	exit;
#}

# set flags as default
my (%flags) = (
"file_input"	=> "false",
"file_output" => "false",
"codon_display"	=> "false"
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
while( defined($temp = shift) ) {
    
    # input file
    if($temp =~ /^[^-]\S*/){
    	$input_file = $temp;
    	$flags{"file_input"} = "true";
    }
    
	# -o : output file
	if($temp =~ /^-\S*o\S*/){
		$flags{"file_output"} = "true";
		$output_file = shift;
		open(OUTPUT_FILE, ">", $output_file) || die "Cannot open:\n\t$!";
	}

	# -t[Number] : get specified number of trials
	if($temp =~ /^-\S*t(\d*)\S*/){
		$numOfTranslation = $1;
		unless (defined $numOfTranslation) { $numOfTranslation = 5; }
	}

	# -c : display codons (dna sequence)
	if($temp =~ /^-\S*c\S*/){
		$flags{"codon_display"} = "true";
	}
	
	# Wrong options
	if($temp =~ /^-[cto]*([^cto\d])+[cto]*$/){
		print "Wrong option used: $1\n";
		usage();
		exit;
	}
}

################################################### PROCESS START #############
print "\n";

####################### Read Input file
if($flags{"file_input"} eq "true") {
	
	PrintOrDisp ("Loading \"$input_file\" ... ");

	open(INPUT_FILE, $input_file) || die "Cannot open: $input_file\n\t$!";
	my $seq = '';
	while( <INPUT_FILE> ) {
		chomp;
		if ($_ =~ /^[#^>]/) {
			# do nothing for the comment line.
		}else{
			$seq = $seq . $_;
		}
	}

	# processing dna sequence for calculation.
	$seq =~ s/\s//g;					# Remove spaces
	$seq =~ s/[^A-Za-z]//g;				# Remove non-alphabet
	$seq =~ tr/[a-z]/[A-Z]/;			# To Upper Character
	#$seq =~ s/[^AGCTNYRSWHVDB]//ig;	# Remove everything except DNA

	close(INPUT_FILE);

	# set input sequence to default multi-codon sequence
	$MultiCodonSequence = $seq;
	$numOfNucleotides = length $MultiCodonSequence;
	$numOfCodons = scalar ($numOfNucleotides / 3);

	#PrintOrDisp ( "$numOfNucleotides of Nucleotides ($numOfCodons codons) were read.\n");
	PrintOrDisp (" Done.\n");
}

####################### Calc.
PrintOrDisp ("Processing ... ");

my @possibleDNASequences = ();		# depository of DNA sequence
my @possiblePeptideSequences = ();	# depository of Peptide sequence

for (my $repeat=0; $repeat<$numOfTranslation; $repeat++) {

	# generate DNA sequence
	my $possibleDNASequence = getPossibleSequence($MultiCodonSequence);
	push(@possibleDNASequences, $possibleDNASequence);

	# generate Peptide sequence
	my $possiblePeptideSequence = dna2peptide($possibleDNASequence);
	push(@possiblePeptideSequences, $possiblePeptideSequence);
}

PrintOrDisp (" Done.\n");

#################################### Report
PrintOrDisp( "\n<Report>\n" );
PrintOrDisp( "=====================================================================\n" );

# display input dna sequence
my $seqByThreeLetters = GetSeqByThreeLetters($MultiCodonSequence);
PrintOrDisp ("[Input Multi-codon DNA Sequence] $numOfNucleotides bases ($numOfCodons codons)\n");
PrintOrDisp ("$seqByThreeLetters\n");


# display generated dna sequence from multi-codons
if ($flags{'codon_display'} eq 'true') {
	
	PrintOrDisp ("\n[DNA Sequences used for each translation]");
	for(my $repeat=0; $repeat<@possibleDNASequences; $repeat++) {

		my $seqByThreeLetters = GetSeqByThreeLetters($possibleDNASequences[$repeat]);
		
		PrintOrDisp ("\n");
		PrintOrDisp(">DNA_Seq%d\n", $repeat);
		PrintOrDisp("$seqByThreeLetters\n");
	}
}

# display translated peptide
PrintOrDisp ("\n[Peptide Sequences]");
for(my $repeat=0; $repeat<@possiblePeptideSequences; $repeat++) {
	PrintOrDisp ("\n");
	PrintOrDisp(">Pep%d\n", $repeat);
	PrintOrDisp("$possiblePeptideSequences[$repeat]\n");
}
print "\n";

# finish program
if ( $flags{'file_output'} eq 'true' ) {
	close OUTPUT_FILE;
}

##################################### End of Main #####################################



##################################### Sub functions ############################
####################### SUB TestTranslation
sub getPossibleSequence {
	my $sampleDNASeq = $_[0];	# argument given as multi-codon DNA sequence

	# creat a DNA sequence from given multi-codon sequence
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

	return $sampleDNASeq;
}

######################## Generate all codons from a single multi-codon
# sub get codons from given a single multicodon
sub getCodonsFromMulticodon {
    my ($aMulticodon) = @_;
    my @codons = ();       # save all codons here
    my $aCodon = '';       # for calc.
    
    # get multibases on each position of codon
    # 만약 멀티코돈이 정의 되지 않았거나 길이가 짧으면 에러를 발생.
    if((!defined $aMulticodon) || (length $aMulticodon) != 3) {
        print "Warning: Given Codon doesn't have 3 nucleotides (or doesn't exist)\n";
    }
    my $thirdBase = substr($aMulticodon, 2, 1);
    my $secondBase = substr($aMulticodon, 1, 1);
    my $firstBase = substr($aMulticodon, 0, 1);

    # get all possible nucleotides at each base
    my @firstBases = getAllPossibleBasesFromMultibase ($firstBase);
    my @secondBases = getAllPossibleBasesFromMultibase ($secondBase);
    my @thirdBases = getAllPossibleBasesFromMultibase ($thirdBase);
    
    # make codons by combinating possible nucleotide sets.
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
sub getAllPossibleBasesFromMultibase {
    my ($aMultibase) = @_;
    
    my @allBases = split(//, $basesFromMultibases{$aMultibase});
    
    #&read_array(\@allBases);
    return @allBases;
}

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
	if ( $flags{'file_output'} eq 'true' ) {
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
    my $someProbability = 0.0;
    
    if( !defined($probability[0]) ){
    	return $nucleotide[ rand(@nucleotide) ];	# calc. using same probability for each nucleotide

    # if probability is given with strong-bond / weak-bond system
    }elsif( scalar(@probability) == 2 ){
		$AT = shift(@probability);
		$CG = shift(@probability);
		unless (($AT + $CG) == 1){
			die "not correct probability\n";
		}

		$someProbability = rand(1);

		if($someProbability < $AT){
			#pick A or T
    		if( rand(1) < 0.5){
    			return 'A';
    		}else{ return 'T'; }

    	}else{
    		# pick C or G
    		if( rand(1) < 0.5 ){
    			return 'C';
    		}else{ return 'G';}
    	}

	# if probability is given for each nucleotide seperately
	}elsif( scalar(@probability) == 4 ){
		($a,$c,$g,$t) = @probability;
		unless( ($a+$c+$g+$t) == 1 ){
			die "not correct probability\n";
		}

		$someProbability = rand(1);

		if($someProbability < $a){
    		return 'A';}
    	elsif($someProbability < $a+$c){
    		return 'C';}
    	elsif($someProbability < $a+$c+$g){
    		return 'G';}
    	else{
    		return 'T';
    	}

	# not allowed probability set was given
	}else{
		die "not correct probability\n";
	}
}

####################### SUB rand_nucleotide_Xs
# This sub routine generates possiblity for each nucleotide and pass that into rand_nucleotide(), finally gets nucleotide according to generated probability set.
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
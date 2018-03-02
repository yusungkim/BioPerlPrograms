####### BackTranslator.pl
####### 16.July.2012 YuSung Kim
####### usage: BackTranslator.pl [input_file]
####### 하는일 : 주어진 단백질 배열을 가지고, DNA배열을 작성한다.
#######
####### 알고리즘 : 1. codon usage의 확률에 따라 아미노산을 DNA로 변경한다.
#######          2. Enzyme site를 검사하고 제거한다. Enzyme file이 필요
#######          3. 아직 실장하지 않았지만, GC content, CAI 등.고려
#######			 4. 5'부분의 2차구조 파괴
#######			 5. 고른 AC/GC　분포


sub version {
    print "\n";
    print "#############################################################\n";
    print "##                Back Translator 1.0                      ##\n";
    print "##     .............................................       ##\n";
    print "##             All Copywrights are reserved by             ##\n";
    print "##           Yusung KIM, 2012.7.16 ~ 2013.04.11             ##\n";
    print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com    ##\n";
    print "#############################################################\n";
}
sub usage{
    print "\nusage: perl $0 <input file> [options] [-o Output_File] [-e Enzyme]\n";
    print "ex)\n perl $0 INPUT.txt -o result.txt -e KpnI -e HindIII\n";
    print "options:\n";
    print "	-v	-> display the process\n";
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

srand();

#my $input_AA = 'PIGDGPVLLP'; #GFP loop
#my $input_AA = 'LSKDPNEKRD'; #GFP loop
#my $input_AA = 'NADFLLRVSGMSMKDIG'; #LexA C'-part of ClpX recognition region
my $peptideInquiry = '';
my $numOfAA = 0;

## usages in 1000 codons
my (@Usages_STOP) = (0.3, 1.0, 2.0);
my (@Codons_STOP) = ('TAG', 'TGA', 'TAA');
my (@Usages_STOPm) = (2.0);
my (@Codons_STOPm) = ('TAA');
my (@Usages_R) = (2.1, 3.6, 3.8, 5.9, 19.7, 20.0);
my (@Codons_R) = ('AGG', 'AGA', 'CGA', 'CGG', 'CGC', 'CGT');
my (@Usages_Rrr) = (19.7, 20.0);	# are codon
my (@Codons_Rrr) = ('CGC', 'CGT');
my (@Usages_Rm) = (10.0, 20.0);	# GC content and rare codon
my (@Codons_Rm) = ('CGC', 'CGT');
my (@Usages_L) = (4.2, 10.2, 11.9, 13.0, 14.3, 48.4);
my (@Codons_L) = ('CTA', 'CTC', 'CTT', 'TTG', 'TTA', 'CTG');
my (@Usages_Lrr) = (10.2, 11.9, 13.0, 14.3, 48.4);
my (@Codons_Lrr) = ('CTC', 'CTT', 'TTG', 'TTA', 'CTG');
my (@Usages_Lm) = (15, 15, 40);
my (@Codons_Lm) = ('TTG', 'TTA', 'CTG');
my (@Usages_C) = (5.2, 6.1);
my (@Codons_C) = ('TGT', 'TGC');
my (@Usages_P) = (5.4, 7.5, 8.6, 20.9);
my (@Codons_P) = ('CCC', 'CCT', 'CCA', 'CCG');
my (@Usages_Pm) = (20.9);
my (@Codons_Pm) = ('CCG');
my (@Usages_I) = (6.8, 23.7, 29.8);
my (@Codons_I) = ('ATA', 'ATC', 'ATT');
my (@Usages_Im) = (23.7, 29.8);
my (@Codons_Im) = ('ATC', 'ATT');
my (@Usages_S) = (8.5, 8.9, 9.1, 9.9, 10.4, 15.2);
my (@Codons_S) = ('TCG', 'TCA', 'TCC', 'AGT', 'TCT', 'AGC');
my (@Usages_Sm) = (10, 10, 20, 20);
my (@Codons_Sm) = ('TCC', 'AGT', 'TCT', 'AGC');
my (@Usages_T) = (9.3, 10.3, 13.7, 22.0);
my (@Codons_T) = ('ACA', 'ACT','ACG', 'ACC');
my (@Usages_Tm) = (5, 13.7, 22.0);
my (@Codons_Tm) = ('ACT','ACG', 'ACC');
my (@Usages_H) = (9.3, 12.5);
my (@Codons_H) = ('CAC', 'CAT');
my (@Usages_Hm) = (10, 20);
my (@Codons_Hm) = ('CAC', 'CAT');
my (@Usages_G) = (9.5, 11.3, 25.5, 27.1);
my (@Codons_G) = ('GGA', 'GGG', 'GGT', 'GGC');
my (@Usages_Gm) = (35, 15);
my (@Codons_Gm) = ('GGT', 'GGC');
my (@Usages_V) = (11.6, 14.3, 19.8, 24.4);
my (@Codons_V) = ('GTA', 'GTC', 'GTT', 'GTG');
my (@Usages_Vm) = (10, 10, 25, 25);
my (@Codons_Vm) = ('GTA', 'GTC', 'GTT', 'GTG');
my (@Usages_Y) = (12.2, 17.5);
my (@Codons_Y) = ('TAC', 'TAT');
my (@Usages_Ym) = (15, 15);
my (@Codons_Ym) = ('TAC', 'TAT');
my (@Usages_K) = (12.4, 35.3);
my (@Codons_K) = ('AAG', 'AAA');
my (@Usages_Km) = (15, 30);
my (@Codons_Km) = ('AAG', 'AAA');
my (@Usages_W) = (13.9);
my (@Codons_W) = ('TGG');
my (@Usages_Q) = (14.6, 28.4);
my (@Codons_Q) = ('CAA', 'CAG');
my (@Usages_Qm) = (20, 30);
my (@Codons_Qm) = ('CAA', 'CAG');
my (@Usages_F) = (16.0, 22.1);
my (@Codons_F) = ('TTC', 'TTT');
my (@Usages_Fm) = (20, 20);
my (@Codons_Fm) = ('TTC', 'TTT');
my (@Usages_A) = (17.1, 21.2, 24.2, 30.1);
my (@Codons_A) = ('GCT', 'GCA', 'GCC', 'GCG');
my (@Usages_Am) = (20, 20, 10, 10);	# GC content
my (@Codons_Am) = ('GCT', 'GCA', 'GCC', 'GCG');
my (@Usages_E) = (18.7, 39.1);
my (@Codons_E) = ('GAG', 'GAA');
my (@Usages_Em) = (18.7, 39.1);
my (@Codons_Em) = ('GAG', 'GAA');
my (@Usages_D) = (19.2, 32.7);
my (@Codons_D) = ('GAC', 'GAT');
my (@Usages_Dm) = (19.2, 32.7);
my (@Codons_Dm) = ('GAC', 'GAT');
my (@Usages_N) = (20.6, 21.4);
my (@Codons_N) = ('AAT', 'AAC');
my (@Usages_Nm) = (20.6, 21.4);
my (@Codons_Nm) = ('AAT', 'AAC');
my (@Usages_M) = (26.4);
my (@Codons_M) = ('ATG');


my (%enzymeDataBase); #DB read from file, RestructionEnzymeList, ex) "KpnI" => "GGTACC"
my $numOfEnzymeDB;

# investigate these enzyme cut sites by default.
my @enzymes = ("SD_Seq", "NcoI", "BamHI", "NdeI", "EcoRI", "KpnI", "XhoI", "HindIII", "SacI", "BglII");
my @cutPositions;
	
######## DATA

my $input_file = 'INPUT_PeptideSequence.txt';
my $output_file = 'RESULT_BackTranslated.txt'; #default output file
my $enzyme_file = 'RestructionEnzymeList.txt';	#default enzyme file

######## Other parameters


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
"codon_display"	=> "true",
"visual"	=> "false"
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

	# -v : visual option
	if($temp =~ /^-\S*v\S*/){
		$flags{"visual"} = "true";
	}

	# Wrong options
	if($temp =~ /^-[cvto]*([^ctvo\d])+[cvto]*$/){
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
	$peptideInquiry = $seq;
	$numOfAA = length $peptideInquiry;

	#PrintOrDisp ( "$numOfNucleotides of Nucleotides ($numOfCodons codons) were read.\n");
	PrintOrDisp (" Done.\n");
}


####################### read restruction enzyme from file
PrintOrDisp ("Loading \"$enzyme_file\" ... ");
open(RES, $enzyme_file) || die "Cannot open: $enzyme_file\n\t$!";
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

# add SD seq into the EnzymeDB
$enzymeDataBase{'SD_Seq'} = 'AGGAG';

$numOfEnzymeDB = keys %enzymeDataBase;
PrintOrDisp ( "Done. ($numOfEnzymeDB Restruction Enzymes)\n");


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



# back translate!
my $backTranslatedCodons = peptide2dna ($peptideInquiry);


my @allCutPositions;
my $tries = 0;

################## Repeat check
do {
$tries++;

# Initialize
# set 'undef' to ref.s for cut positions
undef @allCutPositions;
for (my $j = 0; $j < @cutPositions; $j++) {
	undef $cutPositions[$j];
}

# check whether there is specified enzyme sites on back translated dna.
for(my $j=0; $j<@enzymes; $j++) {
		my ($enzyme) = $enzymes[$j];
		my ($recognitionSite) = $enzymeDataBase{$enzyme};
		my $sizeOfRecognitionSite = length $recognitionSite;
		my $numOfCuts = 0;
		
		if ( $flags{'visual'} eq 'true' ) { PrintOrDisp ("$enzyme: \t$recognitionSite\n"); }
		
		# KMP 방법을 쓰면 더 빠를 수도 있겠지만, monomer가 4개밖에 없는 DNA배열특성상 효율이 크게 향상되지는 않을 것 같다.
		# 그냥 알기 쉽게 차근차근 비교하는 방법으로 search하자.
		my (@backTranslatedCodonsInArray) = split (//, $backTranslatedCodons);
		my (@recognitionSiteInArray) = split (//, $recognitionSite);
		
		# 배열의 처음부터 끝까지 비교
		for (my $ctr1 = 0; $ctr1 < @backTranslatedCodonsInArray - @recognitionSiteInArray; $ctr1++) {
		
			my $isSame = 'false';
			
			# 지정된 효소의 인식부위와 일대일비교
			for (my $ctr2 = 0; $ctr2 < @recognitionSiteInArray; $ctr2++) {
				my $DNA1 = $backTranslatedCodonsInArray[$ctr1 + $ctr2];
				my $DNA2 = $recognitionSiteInArray[$ctr2];
				
				$isSame = $DNA1 eq $DNA2 ? 'true' : 'false';
				
				if($isSame eq 'false') {
					# No Recognition in this region.
					$ctr2 = $sizeOfRecognitionSite;	# Escape the loop			
				}else{
					# Do nothing just Check the next base
				}
			}
			
			if($isSame eq 'true'){ # if recognition site was found)

				# record cut positions at two-dim-array[Index_Enz][Index_Cut_No.] = Recognition site location
				$cutPositions[$j][$numOfCuts] = $ctr1;
				
				# increase cut count
				$numOfCuts++;
			} # end of if
			
		} # end of for()
} # end of searching each Enzyme


# show back translated dna seq.
PrintOrDisp ("\nBack translated ..\n");
PrintOrDisp ("[TryNo.%s] ", $tries);
my $seqByThreeLetters = GetSeqByThreeLetters($backTranslatedCodons);
PrintOrDisp ("$seqByThreeLetters\n");

# visualize enzyme cut site	
PrintOrDisp ("\nChecking unappropriated sequence..\n");
	
for (my $j = 0; $j < @enzymes; $j++) {
		my $cutScreen = '';		# this long string will have and show the cut positions
		my $nextPosition = length $backTranslatedCodons;
		my $indexNext = 0;		# index for @totalCutPositions
	
		# set first cutPosition
		if( defined $cutPositions[$j][0]) {
			$nextPosition = $cutPositions[$j][0];
		}else{
			
			my $enzymeName = $enzymes[$j];
			my $enzymeNameForPrint = "[" . $enzymeName . "]";
			my $recognitionSequence = $enzymeDataBase{$enzymeName};
			
			if ($enzymeName eq "SD_Seq") {
				PrintOrDisp ("%9s (%s)-> None.", $enzymeNameForPrint, $recognitionSequence);
			}else{
				PrintOrDisp ("%9s (%s)-> None.", $enzymeNameForPrint, $recognitionSequence);
			}
			PrintOrDisp("\n");
			#PrintOrDisp (" ( %s )\n", $enzymeDataBase{$enzymes[$j]});
			next;
		}
				
		# print every cut position on the top of seq
		for (my $k = 0; $k< length $backTranslatedCodons; $k++) {
			
			# 인식부위에 도달할 때까지 "-"으로 표시한다.
			if($k<$nextPosition) {
				$cutScreen .= "-";

			# 인식부위에 도달했다.
			}elsif($k == $nextPosition){
				my $recognitionSequence = $enzymeDataBase{$enzymes[$j]};
				
				# 표시할 인식부위(cut location)가 더 남아 있을 경우
				if(++$indexNext < $#{$cutPositions[$j]}+1) {
					$nextPosition = $cutPositions[$j][$indexNext];
					if($k <= $nextPosition - length $recognitionSequence ) {
						$cutScreen .= $recognitionSequence;
						$k = $k - 1 + length $recognitionSequence;
					}else{
						$cutScreen .= "|";	# if another site exist closly behind
					}

				# 더 이상 표시할 인식부위가 없을 경우,
				# 이번 인식부위만 표시해 준 후, 배열의 나머지 부분은 끝까지 "-"로 표시하게끔 설정해 준다.
				}else{

					# 마지막 인식부위를 표시해 준다.
					$cutScreen .= $recognitionSequence;

					# 나머지 부분을 "-"로 채우게끔 설정해 준다.
					$nextPosition = length $backTranslatedCodons;
					$k = $k - 1 + length $recognitionSequence;
				}
			}
		}
		
		# print cut position by three letter spaced way
		my $seqByThreeLetters = GetSeqByThreeLetters($cutScreen);
		my $enzymeNameForPrint = "[" . $enzymes[$j] . "]";
		PrintOrDisp ("%9s ", $enzymeNameForPrint);
		PrintOrDisp ("$seqByThreeLetters\n");
}	

# gather all cut sites
for (my $j = 0; $j < @enzymes; $j++) {
	if( defined $cutPositions[$j][0]) {
		for( my $y=0; $y < $#{$cutPositions[$j]} +1; $y++) {
			push (@allCutPositions, $cutPositions[$j][$y]);
		}
	}
}
@allCutPositions = sort @allCutPositions;

# re-back translate for all cut postions
if (@allCutPositions != 0) {
	print "Modifying enzyme istes ... ";
	for (my $i=0; $i<@allCutPositions; $i++) {
		my $aa = '';
		my $dna = '';
		
		my $pos = $allCutPositions[$i];
		
		if ($pos % 3 == 0) { # change 2 aa
			#change amino acid at pos/3+1
			$aa = substr($peptideInquiry, $pos/3, 2);
			$dna = peptide2dna($aa);
			
			#exchange dna seq
			substr($backTranslatedCodons, $pos, 2*3) = $dna;
			
		}else{	# change 3 aa
			$aa = substr($peptideInquiry, $pos/3, 3);
			$dna = peptide2dna($aa);
		
			#exchange dna seq
			my $offset = $pos % 3;
			substr($backTranslatedCodons, $pos-$offset, 3*3) = $dna;
		}
		my $numOfFixation = $i+1;
		print "$numOfFixation Fixed. ";
	}
	print "All Fixed.\n";
}
	
} while (@allCutPositions != 0);

# Check the back translation was done correctly.
# and print the result.
$peptideInquiry =~ s/\*/_/;	# for comparison original with re-translated aa, the stop signal was changed from '*' to '_'
my $retranslated_AA = dna2peptide($backTranslatedCodons);
if($peptideInquiry eq $retranslated_AA) {
	PrintOrDisp ("Back translation was done.\n");
	$backTranslatedCodons = GetSeqByThreeLetters($backTranslatedCodons);
	PrintOrDisp ("$backTranslatedCodons\n");
}else{
	print STDERR "Program error! retranslated AA is differ from original given aa";
	print "INPUT AA:\n\t$peptideInquiry\n";
	print "ReTranslated AA:\n\t$retranslated_AA\n";
	print "Back translated codons:\n\t$backTranslatedCodons\n";
}

exit;




##################################### SUBs #######################################
### sub amino acid -> codon
sub peptide2dna {
	my($peptide) = @_;
    
    # Initialize variables
    my $dna = '';
    
    # Translate each three-base codon to an amino acid, and append to a protein 
    for(my $i=0; $i < length($peptide) ; $i += 1) {
        $dna .= aa2codon( substr($peptide,$i,1) );
    }
    
    return $dna;
}

sub aa2codon {
	my ($aa) = @_;
	$aa = uc $aa;
	
	my @usages = '';
	my @codons = '';
	my @accumUsage = '';
	
	# set usage
	my (%allUsages) = (
		'_' => \@Usages_STOP,
		'*' => \@Usages_STOP,
		'R' => \@Usages_Rm,
		'L' => \@Usages_Lm,
		'C' => \@Usages_C,
		'P' => \@Usages_Pm,
		'I' => \@Usages_Im,
		'S' => \@Usages_Sm,
		'T' => \@Usages_Tm,
		'H' => \@Usages_Hm,
		'G' => \@Usages_Gm,
		'V' => \@Usages_Vm,
		'Y' => \@Usages_Ym,
		'K' => \@Usages_Km,
		'W' => \@Usages_W,
		'Q' => \@Usages_Qm,
		'F' => \@Usages_Fm,
		'A' => \@Usages_Am,
		'E' => \@Usages_Em,
		'D' => \@Usages_Dm,
		'N' => \@Usages_Nm,
		'M' => \@Usages_M,
	);	
	if(exists $allUsages{$aa}) {
		@usages = @{$allUsages{$aa}};
	}else{
		print STDERR "Bad amino acid \"$aa\"";
	}

	# set codons
	my (%allCodons) = (
		'-' => \@Codons_STOP,
		'*' => \@Codons_STOP,
		'R' => \@Codons_Rm,
		'L' => \@Codons_Lm,
		'C' => \@Codons_C,
		'P' => \@Codons_Pm,
		'I' => \@Codons_Im,
		'S' => \@Codons_Sm,
		'T' => \@Codons_Tm,
		'H' => \@Codons_Hm,
		'G' => \@Codons_Gm,
		'V' => \@Codons_Vm,
		'Y' => \@Codons_Ym,
		'K' => \@Codons_Km,
		'W' => \@Codons_W,
		'Q' => \@Codons_Qm,
		'F' => \@Codons_Fm,
		'A' => \@Codons_Am,
		'E' => \@Codons_Em,
		'D' => \@Codons_Dm,
		'N' => \@Codons_Nm,
		'M' => \@Codons_M,
	);
	if(exists $allCodons{$aa}) {
		@codons = @{$allCodons{$aa}};
	}else{
		print STDERR "Bad amino acid \"$aa\"";
	}
	
	#calc.
	my $totalUsage = 0.0;
	my $selectionRef =0;
	for (my $i = 0; $i < @usages; $i++) {
		$totalUsage += $usages[$i];
		$accumUsage[$i] = $totalUsage;
	}
	my $aUsage = rand($totalUsage);
	my $aCodon = '';
	for (my $i = 0; $i < @usages; $i++) {
		if($aUsage <= $accumUsage[$i]) {
			$aCodon = $codons[$i];
			$selectionRef = $i; # for 'visual' option
			$i = @usages;	# end loop
		}
	}
	
	if ($flags{"visual"} eq "true") {
		print "$aa : ";
		print "$codons[$selectionRef]($usages[$selectionRef]) selected from \t";
		&read_array(\@usages);
		print "\n";
	}
	
	return $aCodon;
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
sub PrintOrDisp {
	my (@args) = @_;
	my $format;
	my @arguments;
	my $withArguments= 'false';
	
	if (@args >1 ) {
		$format = shift;
		(@arguments) = (@args);
		$withArguments = 'true';
	}else{
		$format = shift;
	}
	
	#print in file
	if ( $flags{'file_output'} eq 'true' ) {
		if ($withArguments eq 'true') {
			printf OUTPUT_FILE $format, @arguments;
		}else{
			printf OUTPUT_FILE $format;
		}
		
    #display in stdout
	}else{
		if ($withArguments eq 'true') {
			printf $format, @arguments;
		}else{
			printf $format;
		}
	}
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
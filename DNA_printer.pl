####### DNA_printer.pl
####### 2013.09.24 YuSung Kim
####### usage: DNA_printer [input_file]
####### 하는일 : 주어진 배열을 3개씩 (코돈으로 표시한다)

sub version {
    print "\n";
    print "#############################################################\n";
    print "##                DNA Printer     1.0                      ##\n";
    print "##     .............................................       ##\n";
    print "##             All Copywrights are reserved by             ##\n";
    print "##           Yusung KIM, 2013.9.24 ~ 2013.9.24             ##\n";
    print "##  jinliuxing\@gmail.com, yusungkimlovesyou\@hotmail.com    ##\n";
    print "#############################################################\n";
}
sub usage{
    print "\nusage: perl $0 [input_file] [options]\n";
    print "ex)\n perl $0 input.txt -1\n";
    print "options:\n";
    print "	-<digit>	-> display using given offset digit.\n";
    print " -t          -> translate & display peptide.\n";
    print " -c<digit>   -> number of column to display. (use with -t option)\n";
    print " -t2         -> translate & display peptide mode2.\n";
    print "	-h/?        -> display this\n";
}


################################### MAIN #################################

#!/usr/bin/perl -w
use warnings;
use strict;

my $input_file;
my $output_file;

# this seq is NheI-GFP-KpnI
my $input_DNA = "GCTAGCAAAGGTGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCCGTGGAGAGGGTGAAGGTGATGCTACAATCGGAAAACTCACCCTTAAATTTATTTGCACTACTGGAAAACTGCCTGTTCCGTGGCCAACACTTGTCACTACTCTGACCTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCACATGAAACGTCATGACTTTTTCAAGAGTGCCATGCCGGAAGGTTATGTACAGGAACGCACTATTTCTTTCAAAGATGACGGGAAATACAAGACGCGTGCTGTAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAGGGTACTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAATACAACTTTAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCACAGTTCGCCACAACGTTGAAGATGGTTCCGTTCAACTTGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCTACACAAACTGTCCTTTCGAAAGATCCGAACGAAAAGCGGGACCACATGGTGCTGCACGAGTACGTGAACGCCGCTGGAATCACACATGGTGGTACC";

my $codon_offset = 0;

####### CHECKING OPTIONS
# set flags as default
my (%flags) = (
"visual" => "false", ## donot display the process.
"translate" => "false",
"numberOfColumn" => "5",    #default
"output" => "false",
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
    
    # -v : visual option
    if($temp =~ /^-\S*v\S*/){
        $flags{"visual"} = "true";
    }
    
    # -t : translate & display peptide
    if($temp =~ /^-\S*t\S*/){
        $flags{"translate"} = "true";
    }

    # -c<digit> : number of column
    if ($temp =~ /^-\S*c(\d)\S*/){
        $flags{"numberOfColumn"} = $1;
    }

    # -o : output file
    if($temp =~ /^-\S*o\S*/){
        $flags{"output"} = "true";
        $output_file = shift;
        open(OUTPUT_FILE, ">", $output_file) || die "Cannot open:\n\t$!";
    }
    
    # -<digit> : codon offset
    if($temp =~ /^-(\d)\S*/){
        $codon_offset = $1;

        # offset value should be 0, 1, or 3
        $codon_offset = $codon_offset % 3;
    }
    # Wrong options
    if($temp =~ /^-[voe]*([^toe\d])+[voe]*$/){
        print "Wrong option used: $1\n";
        usage();
        exit;
    }
}


# processing the DNA for display
my $seqByThreeLetters = GetSeqByThreeLetters($input_DNA, $codon_offset);


if ($flags{"translate"} eq "true") {
    my $peptide = dna2peptide($input_DNA);
    my $peptideBy10Letters = peptidePrinter($peptide, $flags{"numberOfColumn"});

    # print peptide
    PrintOrDisp($peptideBy10Letters);

    # print peptide + DNA
    my $scale;
    my @peptideDisplayArray;

    # get first line as scale and remove it from the given seq.
    $peptideBy10Letters =~ s/^(.+)\n//;
    $scale = $1;

    $scale =~ s/ //g;
    $scale =~ s/(\w)/$1   /g;
    $scale = "    ".$scale;

    $peptideBy10Letters =~ s/ //g;
    $peptideBy10Letters =~ s/(\d+)/$1 /g;
    $peptideBy10Letters =~ s/([a-zA-Z])/$1   /g;
    #$peptideBy10Letters =~ s/ /                              /g;
    #$peptideBy10Letters =~ s/^(\d)\s+(\d)\s+(\d)\s*/$1$2$3 /g;
    
    @peptideDisplayArray = split ("\n", $peptideBy10Letters);

    my $dnaXCodonsPerLine = GetXCodonsPerLine($seqByThreeLetters, $flags{"numberOfColumn"}*10);
    my @dnaDisplayArray = split ("\n", $dnaXCodonsPerLine);

    my $display = $scale."\n";
    for(my $i=0; $i<@peptideDisplayArray; $i++) {
        $display .= $peptideDisplayArray[$i]."\n    ".$dnaDisplayArray[$i]."\n";
    }
    PrintOrDisp("\n".$display);


}else{
    PrintOrDisp("%s\n", $seqByThreeLetters);    
}




###################################################### SUB Functions
sub GetXCodonsPerLine {
    my $seq = shift;
    my $numOfCodonsPerLine = shift;

    # if $seq is seqByThreeLetters
    if ($seq =~ /\w\w\w \w\w\w/) {
        # do nothing
    }else{
        # make seqByThreeLetters
        $seq = GetSeqByThreeLetters($seq);
    }

    my $length = length $seq;

    my $numberOfLines = ($length + 1) / ((3+1) * $numOfCodonsPerLine);

    my $result = '';

    my $offset = 0;
    for(my $i=1; $i < $numberOfLines; $i++) {
        $result .= substr($seq, $offset, (3+1)*$numOfCodonsPerLine)."\n";
        $offset = ($i) * (3+1) * $numOfCodonsPerLine;
    }

    $result .= substr($seq, $offset, $length-$offset);

    return $result;

}
sub GetSeqByThreeLetters {
    my $seq = shift;
    my $offset = shift;

    my $length = length $seq;
    my $seqByThreeLetters = '';

    if ($length < 1) {
        return;
    }

    # offset
    if ($offset > 0) {
        $seqByThreeLetters = substr($seq, 0, $offset);
        $seqByThreeLetters .= ' ';

        $seq = substr($seq, $offset, $length-$offset);
        $length = length $seq;
    }
    
    # 3 letter dividing
    for(my $i=0; $i<$length-2; $i=$i+3) {
        $seqByThreeLetters .= substr($seq, $i, 3);
        $seqByThreeLetters .= ' ';
    }

    # if last 1 or 2 letters were left for printing
    if ($length % 3 == 1) {
        $seqByThreeLetters .= substr($seq, -1);
    }elsif ($length % 3 == 2) {
        $seqByThreeLetters .= substr($seq, -2);
    }

    chomp $seqByThreeLetters;

    return $seqByThreeLetters;
}

sub peptidePrinter {
    my $peptide = shift;
    my $numberOfColumn = shift;
    my $pepLength = length ($peptide);

    
    #$peptideBy10Letters =  "    1/51       11/61      21/71      31/81      41/91     \n";
    #$peptideBy10Letters .= "    1234567890 1234567890 1234567890 1234567890 1234567890\n";
    #$peptideBy10Letters .= "001 ";

    my $peptideBy10Letters = "   ";
    for(1..$numberOfColumn) {
        $peptideBy10Letters .= " 1234567890";
    }
    $peptideBy10Letters .= "\n";

    my $i=1;
    my $leftCounter = '';

    # write 50 peptide
    for($i=1; $i * $numberOfColumn * 10 <= $pepLength; $i++){

        # write left counter
        if( ($i-1) * $numberOfColumn * 10 < 10) {
            $leftCounter = "001 ";
        }elsif( ($i-1) * $numberOfColumn * 10 < 100) {
            $leftCounter = "0".(($i-1)*$numberOfColumn*10+1)." ";
        }else{
            $leftCounter = ($i-1) *$numberOfColumn *10 +1;
            $leftCounter .= " ";
        }
        $peptideBy10Letters .= $leftCounter;

        # write peptide
        for(my $k=0; $k<$numberOfColumn; $k++) {
            $peptideBy10Letters .= substr($peptide, ($i-1)*$numberOfColumn *10 + $k*10, 10);
            if ($k==$numberOfColumn-1) {
                $peptideBy10Letters .= "\n";
            }else{
                $peptideBy10Letters .= " ";
            }
        }
    }

    # write less than 50 peptides or the rest of the peptide
    if( ($i-1) * $numberOfColumn * 10 < 10) {
        $leftCounter = "001 ";
    }elsif( ($i-1) * $numberOfColumn * 10 < 100) {
        $leftCounter = "0".(($i-1)*$numberOfColumn*10+1)." ";
    }else{
        $leftCounter = ($i-1) *$numberOfColumn *10 +1;
        $leftCounter .= " ";
    }
    $peptideBy10Letters .= $leftCounter;

    for(my $k=0; $pepLength - ($i-1)*$numberOfColumn*10 >= ($k+1)*10; $k++) {
        $peptideBy10Letters .= substr($peptide, ($i-1)*$numberOfColumn*10 + $k*10, 10);

        if ($pepLength - ($i-1)*$numberOfColumn*10 == ($k+1)*10 ) {
            $peptideBy10Letters .= "\n";
        }else{
            $peptideBy10Letters .= " ";
        }
    }

    my $theNumberOfTheRestPeptide = $pepLength % 10;
    $peptideBy10Letters .= substr($peptide, -$theNumberOfTheRestPeptide);
    $peptideBy10Letters .= "\n";

    return $peptideBy10Letters;
}


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
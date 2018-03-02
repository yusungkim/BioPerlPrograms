sub copyright {
  print "####################################################\n";
  print "##    PrimerGenerator Version 1.0                 ##\n";
  print "##  ............................................. ##\n";
  print "##     using Smith-Waterman algorithm. O(n^2)     ##\n";
  print "##        Kim Yu-Sung, 2011.09.14 ~ .             ##\n";
  print "####################################################\n";
}
sub usage{
  print " usage: perl $0 [OPTIONS..] SUBJECT QUERY\n";
  print " options:
  		-v Version\n
		-gNUM GAP_PENALTY\n
		-mSTR MATCH_MATRIX\n
		#Caution : identity is calculated under non-gap alignment.
		-iNUM IDENTITY_SETTING 'NUM' PERSENTAGE 65, 75, 95, 99\%\n";
}

#!/usr/bin/perl -w
use warnings;
use strict;

my ($i, $j);

my ($randomNucleotide) = rand_nucleotide_N();
my ($multimer10) = "GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN GTN CCN GGN GTN GAR GGN GTN CCN GGN GTN GTN GGN";
my ($maxPrimerSize) = 200;
my ($primerSize) = 0;
my ($overLapSize) = 20;
my ($length) = length ($multimer10);

# remove spaces
$multimer10 =~ s/\s//g;

# split into 60 to 200 size of window (primer)

#print "$multimer10\n";
print " > size of multimer (length) = $length\n";

my ($r) = 0;
my ($n) = 0;

my ($flagOut) = 0;
my ($divider) = 0;
for ($i=0; $flagOut == 0;$i++) {
	$divider = ($maxPrimerSize - $overLapSize -$i);
	$r = ($length - $overLapSize) % $divider;
	$n = int($length / $divider);
	
	print " > r is $r \t n is $n \t divider is $divider\n";
	
	if ($r == 0) { $flagOut = 1; }
		
	if ($overLapSize + $r <= $n * $i) { $flagOut = 1; }
}
#print " > r is $r \t n is $n \t divider is $divider\n";

my (@primers);
if ($r == 0) {
	$primerSize = $maxPrimerSize - $i;
	my ($position) = 0;
	for($j = 0; $position + $overLapSize < $length; $j++) {
		$primers[$j] = substr $multimer10, $position, $primerSize;
		$position += ($primerSize - $overLapSize);
	}
}else{
	my ($lengthForDivide) = $length + ($n - 1) * $overLapSize;
	$primerSize = int ($lengthForDivide / $n);
	$r = $lengthForDivide % $n;
	print " > r is $r \t n is $n \t primerSize is $primerSize\n";
	
	my ($position) = 0;
	for($j = 0; $position + $overLapSize < $length; $j++) {
		if ($r == 0) {
			$primers[$j] = substr ($multimer10, $position, $primerSize);
			$position += ($primerSize - $overLapSize);
			
		}else {
			$primers[$j] = substr ($multimer10, $position, $primerSize+1);
			$position += ($primerSize+1 - $overLapSize);
			$r--;
		}
	}
}

print scalar(@primers);
print "\n";
print length(@primers);


#for each primer {
 #for x times {
#	random seq. front and rear
#	save expted score
#	compair with all seq

	


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
# N : A, C, G, T
sub rand_nucleotide_N {
	return rand_nucleotide (0.25, 0.25, 0.25, 0.25);
}
# Y : pYrimidine, C, T
sub rand_nucleotide_Y {
	return rand_nucleotide (0.0, 0.5, 0.0, 0.5);
}
# R : puRine, A, G
sub rand_nucleotide_R {
	return rand_nucleotide (0.5, 0.0, 0.5, 0.0);
}
# S : Strong, C, G
sub rand_nucleotide_S {
	return rand_nucleotide (0.0, 0.0, 0.5, 0.5);
}
# W : Week, A, T
sub rand_nucleotide_W {
	return rand_nucleotide (0.5, 0.5, 0.0, 0.0);
}
# H : A, C, T
sub rand_nucleotide_H {
	return rand_nucleotide (0.33333, 0.33333, 0.0, 0.33334);
}
# V : A, C, G
sub rand_nucleotide_V {
	return rand_nucleotide (0.33333, 0.0, 0.33333, 0.33334);
}
# D : A, G, T
sub rand_nucleotide_D {
	return rand_nucleotide (0.33333, 0.33333, 0.33334, 0.0);
}
# B : C, G, T
sub rand_nucleotide_B {
	return rand_nucleotide (0.0, 0.33333, 0.33333, 0.33334);
}
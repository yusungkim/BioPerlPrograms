#! /usr/bin/perl
### Single Molecular Experiment _ File Concatanator.
### 24. Sep. 2008,  Yu-Sung Kim
### ~ 25. Sep. 2008 Version 1.1

### perl SME_concatanator.pl <file-1.txt>

## <file-1> should be "xxxx-1.txt" or "xyz-1.txt"
## then all sirese of files, "xxxx-1.txt", "xxxx-2.txt", ... will be concatenated.

## inside of "file-1.txt" should be like this...
## <headers>:non-disit letter start
## Scene	X	Y	X	Y	Time	Angle	Revolution	Error
## 1	12.9349	23.1233	0.0000	30.0923 4.9938	0
## 2 12.9349	23.1233	0.0000	30.0923 4.9938	0
## <end of file>

use warnings;
use strict;

my $fileName_INPUT1 = '';
my $fileName_INPUT2 = '';
my $fileName_OUTPUT = ''; ## added fille
my $scene;
my $scene2;
my $x_ROI;
my $y_ROI;
my $x_ALL;
my $y_ALL;
my $time;
my $time2;
my $angle;
my $revolution;
my $revolution2;
my $error;

my $line2;
my $dummy;

## get filename from command line
if(@ARGV < 1){ usage(); exit; }
$fileName_INPUT1 = <@ARGV>;

## making OUTPUT file name
if ($fileName_INPUT1 =~ /^(.*)\d.txt$/) {
	$fileName_OUTPUT = $1."all.txt";
} else {
	die "INPUT file_name type error!\n";
}


### READ FILE DATA
open(INPUT1, "$fileName_INPUT1") || die "Cannot open:\n\t$!";
open(OUTPUT, "+> $fileName_OUTPUT") || die "Cannot open:\n\t$!";
print "Creating concatanated file ...";
print "\n  ... Adding $fileName_INPUT1";

## discard headers
while ( !((my $line =<INPUT1>) =~ /^Scene|^Frame/i) ) { ; }

## write output
print OUTPUT "Scene\tX_ROI\tY_ROI\tX_ALL\tY_ALL\ttime\tangle\trevolution\r\n";
while ( my $line = <INPUT1> ) {
		$line =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
		$scene = $1;
		$x_ROI = $2;
		$y_ROI = $3;
		$x_ALL = $4;
		$y_ALL = $5;
		$time = $6;
		$angle = $7;
		$revolution = $8;
		#$error = $9;
		
		printf OUTPUT ("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\r\n",$scene,$x_ROI,$y_ROI, $x_ALL, $y_ALL, $time, $angle, $revolution);
		#print OUTPUT "$scene\t$x_ROI\t$y_ROI\t$x_ALL\t$y_ALL\t$time\t$angle\t$revolution\r";
}
close(INPUT1) || die "Cannot close: \n\t$!";
print " ... Done.";

### Concatenate files
## create next file name
do {
	my $next;
	$fileName_INPUT1 =~ /^(.*-)(\d+).txt$/;
	$fileName_INPUT2 = $1;
	$next = $2 + 1;
	$fileName_INPUT2 .= $next.".txt";
	
	if (-e $fileName_INPUT2) {
		## concatanate
		print "\n  ... Adding $fileName_INPUT2";
		open(INPUT2, $fileName_INPUT2) || die "Cannot open $fileName_INPUT2\n\t$!";
		
		## discard headers
		while ( !((my $line =<INPUT2>) =~ /^Scene|^Frame/i) ) { ; }
		
		## adding data
		while( <INPUT2> ) {
			$_ =~ /^([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+).*/;
			($scene2, $x_ROI, $y_ROI, $x_ALL, $y_ALL, $time2, $angle, $revolution2) = ($1, $2, $3, $4, $5, $6, $7, $8);
			
			## shift data to make gradual increasement
			$scene2 += $scene;	## $scene has last data of previous file.
			$time2 += $time;
			$revolution2 += $revolution;
			#print OUTPUT "$scene2\t$x_ROI\t$y_ROI\t$x_ALL\t$y_ALL\t$time2\t$angle\t$revolution2\r\n";	
			printf OUTPUT ("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\r\n",$scene2,$x_ROI,$y_ROI, $x_ALL, $y_ALL, $time2, $angle, $revolution2);
		} ## end of while
		close(INPUT2) || die "Cannot close $fileName_INPUT2\n\t$!";
		print " ... Done.";
		
	} ## end of if
	
	## prepare for the next roop
	$fileName_INPUT1 = $fileName_INPUT2;
	$scene = $scene2;	## save last data
	$time = $time2;
	$revolution = $revolution2;
	
} while (-e $fileName_INPUT2);
print "\n";

close(OUTPUT) || die "Cannot open:\n\t$!";
print "  ... $fileName_OUTPUT ... Created\n";

### USAGE
sub usage{
	print "\n";
	print "  Usage  : <file1> <file2\n";
	print "\n";
}
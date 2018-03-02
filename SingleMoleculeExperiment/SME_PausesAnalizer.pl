#! /usr/bin/perl
### Single Molecular Experiment _ Pause Analizer.
### 22. Sep. 2008,  Yu-Sung Kim
### ~ 29. Sep. 2008 Version 1.2

### perl SME_PauseAnalizer.pl <data_file>

## "data_file" should be like this...
## <headers>
## --------------
## Scene	X	Y	X	Y	Time	Angle	Revolution	Error
## 1	12.9349	23.1233	0.0000	30.0923 4.9938	0
## 2 12.9349	23.1233	0.0000	30.0923 4.9938	0
## <end of file>

use warnings;
use strict;

my $separator_frame = 30;
my $separator_rot = 2;
my $rotORstop;
my $fileName = '';
my $processDataFileName1 = '';
my $processDataFileName2 = ''; ## for Revolution Buffer
my $processDataFileName3 = ''; ## for Delta_Rev.
my $processDataFileName4 = ''; ## for Rotation
my $processDataFileName5 = ''; ## for Stop
my $scene;
my $scene2;
my $x_ROI;
my $y_ROI;
my $x_ALL;
my $y_ALL;
my $time;
my $angle;
my $revolution;
my $revolution2;
my $delta_rev;
my $error;

my $line2;
my $dummy;

### OPEN BULK_DATA_File

## get filename from command line
if(@ARGV == 1) {
	$fileName = <@ARGV>;
}elsif(@ARGV == 3){
	($fileName, $separator_rot, $separator_frame) = @ARGV;
}else{
	usage();
	exit;
}

print "Dividing Rotating/Pausing data ... $separator_rot rot / $separator_frame frame";

## making process file names
if ($fileName =~ /^(.*).txt$/) {
	$processDataFileName1 = $processDataFileName2 = $processDataFileName3 = $1;
	$processDataFileName4 = $processDataFileName5 = $1;
} else {
	$processDataFileName1 = $processDataFileName2 = $processDataFileName3 = $fileName;
	$processDataFileName4 = $processDataFileName5 = $fileName;
}
$processDataFileName1 .= ".dummy1.txt";
$processDataFileName2 .= ".dummy2.txt";
$processDataFileName3 .= ".dummy3.txt";
$processDataFileName4 .= ".pause.txt";
$processDataFileName5 .= ".pause_period.txt";


### READ FILE DATA
open(BULKDATA, "$fileName") || die "Cannot open:\n\t$!";

## discard headers
while ( !((my $line =<BULKDATA>) =~ /^(Scene|^Frame).*$/) ) { ; }

### make processDataFile
## Create process Data_file_Name
open(PROCESSDATA1, "+> $processDataFileName1")  || die "Cannot open file:\n\t$!";
print PROCESSDATA1 "Scene\tRotation\n";

open(PROCESSDATA2, "+> $processDataFileName2")  || die "Cannot open file:\n\t$!";
print PROCESSDATA2 "Scene\tTime\tX\tY\tAngle\tRotation\n";

## extract data
foreach my $line (<BULKDATA>) {
	$line =~ /^([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
	$scene = $1;
	$x_ROI = $2;
	$y_ROI = $3;
	$x_ALL = $4;
	$y_ALL = $5;
	$time = $6;
	$angle = $7;
	$revolution = $8;
	#$error = $9;
	print PROCESSDATA1 "$scene\t$revolution\r\n";
	print PROCESSDATA2 "$scene\t$time\t$x_ALL\t$y_ALL\t$angle\t$revolution\r\n";
}
close(PROCESSDATA1) || die "Cannot close:\n\t$!";
print "\n  ... $processDataFileName1 ... Created";
close(PROCESSDATA2) || die "Cannot close:\n\t$!";
print "\n  ... $processDataFileName2 ... Created";
close(BULKDATA) || die "Cannot close:\n\t$!";

### Process 1
open(PROCESSDATA1, "$processDataFileName1")  || die "Cannot open file:\n\t$!";
open(PROCESSDATA2, "$processDataFileName2")  || die "Cannot open file:\n\t$!";
open(PROCESSDATA3, "+> $processDataFileName3")  || die "Cannot open file:\n\t$!";
$dummy = <PROCESSDATA1>; ## skip header
$dummy = <PROCESSDATA2>; ## skip header
print PROCESSDATA3 "Scene\tTime\tX\tY\tAngle\tRotation\tRotOrStop\tDelta_Rot\r\n"; ## add header

## Skip some data
for (my $i = 0; $i < $separator_frame; $i++) {
	my $dummy = <PROCESSDATA1>;
}

foreach my $line (<PROCESSDATA2>) {
	$line =~ /^([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*$/;
	$scene = $1;
	$time = $2;
	$x_ALL = $3;
	$y_ALL = $4;
	$angle = $5;
	$revolution = $6;
	$line2 = <PROCESSDATA1>;
	#print "1\n";
	#if ($scene > 40) {next;}
	if( defined $line2 ) {
		#print "\2\n";
		$line2 =~ /^([-\d.]+)\s+([-\d.]+).*/;
		$scene2 = $1;
		$revolution2 = $2;
		$delta_rev = $revolution2 -$revolution;
		#print "$revolution2 \n";
		if($delta_rev > $separator_rot) {
			$rotORstop = 1;
		}else{
			$rotORstop = 0;
		}
		print PROCESSDATA3 "$scene\t$time\t$x_ALL\t$y_ALL\t$angle\t$revolution\t$rotORstop\t$delta_rev\r\n";
	}
}
close(PROCESSDATA1) || die "Cannot close:\n\t$!";
close(PROCESSDATA2) || die "Cannot close:\n\t$!";
close(PROCESSDATA3) || die "Cannot close:\n\t$!";
print "\n  ... $processDataFileName3 ... Created";

### Process 2
open(PROCESSDATA3, "$processDataFileName3")  || die "Cannot open file:\n\t$!";
open(PROCESSDATA4, "+> $processDataFileName4")  || die "Cannot open file:\n\t$!";

$dummy = <PROCESSDATA3>; ## skip header
print PROCESSDATA4 "Scene\tTime\tStop_X\tStop_Y\tStop_Angle\tStop_Revolution\tRot_X\tRot_Y\tRot_Angle\tRot_Revolution\tRotOrStop\tDelta_Rot\r\n"; ## add header

my $time_of_start = 0;
my $period = 0;
my @period_rot = ();
my @period_stop = ();
my @test_change = ();
my $previous = 0;
my $first_state;

## read first line.
$dummy = <PROCESSDATA3>;
$dummy =~  /^([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*/;
($scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev) = ($1, $2, $3, $4, $5, $6, $7, $8);

## set first state.
$time_of_start = $time;
$first_state = $previous = $rotORstop;

## write first line.
if ($rotORstop eq "0") { ## if stoped
	printf PROCESSDATA4 ("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t\t\t\t\t%d\t%.4f\r\n", $scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev);
}else{
	printf PROCESSDATA4 ("%d\t%.4f\t\t\t\t\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\r\n", $scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev);
}

## process from the 2nd line.
foreach my $line (<PROCESSDATA3>) {
	$line =~  /^([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*/;
	($scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev) = ($1, $2, $3, $4, $5, $6, $7, $8);
	
	if ($previous eq $rotORstop) {
		;
	}else {
		$period = $time - $time_of_start;
		push (@test_change, $scene);
		
		if($rotORstop eq "0") { ## now stopped
			push (@period_rot, $period);
		}else {
			push (@period_stop, $period);
		}
		
		$time_of_start = $time; ## set for next
		$previous = $rotORstop; ## set for next
		#print "$scene\t$period\n";
	}

	if ($rotORstop eq "0") { ## if stoped
		printf PROCESSDATA4 ("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t\t\t\t\t%d\t%.4f\r\n", $scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev);
	}else{
		printf PROCESSDATA4 ("%d\t%.4f\t\t\t\t\t%.4f\t%.4f\t%.4f\t%.4f\t%d\t%.4f\r\n", $scene, $time, $x_ALL, $y_ALL, $angle, $revolution, $rotORstop, $delta_rev);
	}
} ## end of foreach

## remove first data
if ($first_state eq "1" or $first_state == 1) { ## if first state is rotating
	shift @period_rot;	## remove non-profit data (first data)
}else{
	shift @period_stop;
}

close(PROCESSDATA3) || die "Cannot close:\n\t$!";
close(PROCESSDATA4) || die "Cannot close:\n\t$!";
print "\n  ... $processDataFileName4 ... Created";

#print (@period_rot);
#foreach (@test_change) {
#	print "$_\n";
#}
open(PROCESSDATA5, "+> $processDataFileName5")  || die "Cannot open file:\n\t$!";
my $length_stop = @period_stop;
my $length_rot = @period_rot;

unless (defined $length_stop) { $length_stop = 0; }
unless (defined $length_rot) { $length_rot = 0; }

my $length_larger;
my $length_smaller;
($length_stop < $length_rot) ? ($length_larger = $length_rot, $length_smaller = $length_stop) : ($length_larger = $length_stop, $length_smaller = $length_rot);
#print "\nstop = $length_stop\trot = $length_rot\tlarger = $length_larger\tsmaller = $length_smaller\n";

print PROCESSDATA5 "Stop_Period\tRot_Period\r\n";
for (my $i=0; $i<$length_smaller; $i++) {
	printf PROCESSDATA5 ("%.4f\t%.4f\r\n", $period_stop[$i], $period_rot[$i]);
}
for (my $i=$length_smaller; $i<$length_larger; $i++) {
	if($length_stop > $length_rot) {
		printf PROCESSDATA5 ("%.4f\t\r\n", $period_stop[$i]);
	}else {
		printf PROCESSDATA5 ("\t%.4f\r\n", $period_rot[$i]);
	}
}

close(PROCESSDATA5) || die "Cannot close:\n\t$!";
print "\n  ... $processDataFileName5 ... Created\n";

unlink $processDataFileName1;
print "  ... $processDataFileName1 ... Removed\n";
unlink $processDataFileName2;
print "  ... $processDataFileName2 ... Removed\n";
unlink $processDataFileName3;
print "  ... $processDataFileName3 ... Removed\n";

### USAGE
sub usage{
	print "\n";
	print "  Usage  : <file name> <rotation> <frame>\n";
	print "           smaller than <rotation> rot / <frame> frame will be counted as long-pause\n";
	print "\n";
}
#! /usr/bin/perl
### Single Molecular Experiment _ StepViewer_input file generator from ImageJ_output file.
### 3. Oct. 2008,  Yu-Sung Kim
### ~ 3. Oct. 2008 Version 1.1

### perl SME_ImageJ2StepViewer.pl <file_IJ.txt>

## ImageJ output file type
## Frame\tX	Pixcel\tY Pixcel\tAngle\tAccum\tRevo\n
## 1	12.933449	23.144233	0.000000	00.000000 0.000000
## <end of file>

## StepViewer input file type
## Frame\tX\tY\r\n
## < each lane should be same width. >

use warnings;
use strict;

my $usage = "perl SME_ImageJ2StepViewer.pl <file.IJ.txt>\n \tfile list is also available. \"| *IJ.txt\"\n";

if(@ARGV < 1){ print $usage; exit; }

## process all file in the list
foreach (<@ARGV>) {	doit($_); }

sub doit {
	my ($file_Input) = @_; ## ImageJ output file
	my $file_Output = ''; ## StepViewer input file
	my $frame;
	my $x_ROI;
	my $y_ROI;
	my $angle;
	my $accum;
	my $revolution;
	
	## making OUTPUT file name
	if ($file_Input =~ /^(.*)IJ.txt$/) {
		$file_Output = $1."forSV.txt";
	} else {
		die "INPUT file_name type error!\n";
	}


	### READ FILE DATA
	open(INPUT, "$file_Input") || die "Cannot open:\n\t$!";
	open(OUTPUT, "+> $file_Output") || die "Cannot open:\n\t$!";
	print " Creating \"$file_Output\"  ...";

	## discard headers
	while ( !((my $line =<INPUT>) =~ /^Scene|^Frame|^#frame/i) ) { ; }

	## write output
	print OUTPUT "Frame\tX_ROI\tY_ROI\tAngle\Revolution\r\n";
	while ( my $line = <INPUT> ) {
		$line =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
		$frame = $1;
		$x_ROI = $2;
		$y_ROI = $3;
		$angle = $4;
		$accum = $5;
		$revolution = $6;
		
		printf OUTPUT ("%5d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\r\n",$frame,$x_ROI,$y_ROI,$angle,$revolution);
	}
	close(INPUT) || die "Cannot close: \n\t$!";
	close(OUTPUT) || die "Cannot close: \n\t$!";
	print " ... Done.\n";
}
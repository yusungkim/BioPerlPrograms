#! /usr/bin/perl
### Single Molecular Experiment _ StepViewer_input file generator from ImageJ_output file.
### it also corrects drift if there is refference file.
### 3. Oct. 2008,  Yu-Sung Kim
### ~ 3. Oct. 2008 Version 1.1

### perl SME_ImageJ2StepViewer.pl <list of *IJ.txt>

## ImageJ output file type
## Frame\tX	Pixcel\tY Pixcel\tAngle\tAccum\tRevo\n
## 1	12.933449	23.144233	0.000000	00.000000 0.000000
## <end of file>

## StepViewer input file type
## Frame\tX\tY\r\n
## < each lane should be same width. >

use warnings;
use strict;

my $usage = "perl SME_DriftCorrector.pl <*IJ.txt>\n \tfile list is also available. \"*IJ.txt\"\n";

if(@ARGV < 1){ print $usage; exit; }

## process all file in the list
foreach (<@ARGV>) {	doit($_); }

sub doit {
	my ($file_Input) = @_; ## ImageJ output file
	my $file_Ref = '';	## refference for drift correction.
	my $file_Output = ''; ## Drif corrected file. file type is ImageJ output file.
	my $ref_exist = 'no'; ## yes / no
	my $frame;
	my $x_ROI;
	my $y_ROI;
	my $angle;
	my $accum;
	my $revolution;
	
	## making OUTPUT file name
	if ($file_Input =~ /^(.*)IJ.txt$/) {
		$file_Output = $1."forSV.txt";
		$file_Ref = $1."IJ.ref.txt";
	} else {
		die "INPUT file_name type error!\n";
	}


	### OPEN FILEs
	open(INPUT, "$file_Input") || die "Cannot open:\n\t$!";
	if(-e $file_Ref) {
		$ref_exist = "yes";
		open(REF, "$file_Ref") || die "Cannot open:\n\t$!";
	}
	open(OUTPUT, "+> $file_Output") || die "Cannot open:\n\t$!";
	print " Creating \"$file_Output\"  ...";

	## discard header and title
	while ( !((my $line =<INPUT>) =~ /^Scene|^Frame|^#frame/i) ) { ; }
	if( $ref_exist eq 'yes' ) {
		print " Drift Correcting ...";
		while ( !((my $line =<REF>) =~ /^Scene|^Frame|^#frame/i) ) { ; }
	}else {
		print " No refference file ...";
	}
	
	## write output
	
	
	## WHEN Refference file doesn't exist
	if($ref_exist eq 'no') {
		
		## write title
		print OUTPUT "Frame\tX_ROI\tY_ROI\tAngle\tRevolution\r\n";
		
		## wirte data
		while ( my $line = <INPUT> ) {
			$line =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
			($frame, $x_ROI, $y_ROI, $angle, $accum, $revolution) = ($1, $2, $3, $4, $5, $6);
			
			printf OUTPUT ("%5d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\r\n",$frame,$x_ROI,$y_ROI,$angle,$revolution);
		}
	
	## WHEN Refference file does exist
	}else{
		## write title
		print OUTPUT "Frame\tX_ROI\tY_ROI\r\n";
		
		## wirte data
		my $ref;
		while ( my $line = <INPUT>) {
			$line =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
			($frame, $x_ROI, $y_ROI, $angle, $accum, $revolution) = ($1, $2, $3, $4, $5, $6);
			
			$ref = <REF>;
			$ref =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
			$x_ROI -= $2;
			$y_ROI -= $3;
			
			printf OUTPUT ("%5d\t%6.4f\t%6.4f\r\n",$frame,$x_ROI,$y_ROI);
		}
		close(REF) || die "Cannot close: \n\t$!";
	}
	
	close(INPUT) || die "Cannot close: \n\t$!";
	close(OUTPUT) || die "Cannot close: \n\t$!";
	print " Done.\n";
}
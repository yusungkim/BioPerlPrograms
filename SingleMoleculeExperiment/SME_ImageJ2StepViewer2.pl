#! /usr/bin/perl
### Single Molecular Experiment _ StepViewer_input file generator from ImageJ_output file.
### it also corrects drift if there is refference file.
### And Doit against All of sub-directory.
### 3. Oct. 2008,  Yu-Sung Kim
### ~ 8. Oct. 2008 Version 2.1

### Version Info.
### 2.0 + Recursively against all sub directory.
### 2.1 + display corrected.
### 2.2 + display result modified.

### perl SME_ImageJ2StepViewer.pl <list of *IJ.txt> or <directory>

## ImageJ output file type
## Frame\tX	Pixcel\tY Pixcel\tAngle\tAccum\tRevo\n
## 1	12.933449	23.144233	0.000000	00.000000 0.000000
## <end of file>

## StepViewer input file type
## Frame\tX\tY\r\n
## < each lane should be same width. >

use warnings;
use strict;

my $usage = "perl SME_DriftCorrector.pl <*IJ.txt>\n \tfile list or directory is also available. \"*IJ.txt\"\n";

## for display result
my $num_of_dir = 0;
my $num_of_file = 0;
my $num_of_rfd_file = 0;
my $num_of_error = 0;
my @error_files =();

if(@ARGV < 1){ print $usage; exit; }

## process all file in the directory
if (-d $ARGV[0]) { ## if is directory
	doit_directory($ARGV[0]);
	
## process all file in the list
}else {
	foreach (<@ARGV>) {
		doit($_);
		$num_of_file++;
	}
}

## print result overview
if($num_of_dir > 1) {
	print "$num_of_dir directories are surveyed.\n";
}else{
	print "$num_of_dir directory is surveyed.\n";
}

$num_of_rfd_file += $num_of_file;
if($num_of_file > 1) {
	print "$num_of_file \"*.IJ.txt\" files are processed.\n";
}else{
	print "$num_of_file \"*.IJ.txt\" file is processed.\n";
}
if($num_of_file > 1) {
	print "$num_of_rfd_file files are created.\n";
}else{
	print "$num_of_rfd_file file is created.\n";
}
if($num_of_error > 0) {
	print "$num_of_error Error occured.\n";
	foreach my $one_file (@error_files) {
		print "\t$one_file\n";
	}
}else{
	print "No errors. ¥n";
}

######################################### SUB

sub doit_directory {
	my ($dir) = @_;
	my @args =();
	
	$num_of_dir++;
	
	## doit
	open(FILE_LIST, "ls -l $dir |") or die "cannot open:$!";
	print "$dir :\n";
	foreach (<FILE_LIST>) {
		$_ =~ /\s(\S+)\s*$/;
		my $file_name = $dir."/".$1;
		if($file_name =~ /IJ.txt$/ && -f $file_name) {
			doit($file_name);
			$num_of_file++;
		}
	}
	close(FILE_LIST);
	print "\n";
	
	## doit sub_dir recursively
	open(DIR_LIST, "ls -l $dir |") or die "cat not open:$!";
	foreach (<DIR_LIST>) {
		$_ =~ /\s(\S+)\s*$/;
		my $file_name = $dir."/".$1;
		
		## directory
		if(-d $file_name) {
			doit_directory($file_name);
		}
	}
	close(DIR_LIST);
}

sub doit {
	my ($file_Input) = @_; ## ImageJ output file
	my $file_Ref = '';	## refference for drift correction.
	my $file_Output = ''; ## Output file, file type is ImageJ output file.
	my $file_Output_Rfd = ''; ## Output file (Drift corrected)
	my $ref_exist = 'no'; ## yes / no
	my $file_name_only = '';
	my $frame;
	my $x_ROI;
	my $y_ROI;
	my $angle;
	my $accum;
	my $revolution;
	my $ref = '';
	my $counter = 0;
	my $skipped_lines =0;
	
	## making OUTPUT file name
	if ($file_Input =~ /^(.*)IJ.txt$/) {
		$file_Output = $1."forSV.txt";
		$file_Output_Rfd = $1."rfd.forSV.txt";
		$file_Ref = $1."IJ.ref.txt";
	} else {
		die "INPUT file_name type error!\n";
	}


	### OPEN FILEs
	open(INPUT, "$file_Input") || die "Cannot open:\n\t$!";
	open(OUTPUT, "+> $file_Output") || die "Cannot open:\n\t$!";
	
	if(-e $file_Ref) {
		$ref_exist = "yes";
		open(REF, "$file_Ref") || die "Cannot open:\n\t$!";
		open(OUTPUT_RFD, "+>$file_Output_Rfd") || die "Cannot open:\n\t$!";
	}
	
	## for displaying..
	$file_Input =~ /([^\s^\/]+)$/;
	$file_name_only = $1;
	print "\tIN  << \"$file_name_only\"";

	## discard header and title
	while ( !((my $line =<INPUT>) =~ /^Scene|^Frame|^#frame/i) ) { ; }
	if( $ref_exist eq 'yes' ) {
		print " (WITH Refference)\n";
		while ( !((my $line =<REF>) =~ /^Scene|^Frame|^#frame/i) ) { ; }
	}else {
		print " (NO Refference)\n";
	}
	
	## write output	
	## write title
	## Non-Drift corrected file.
	print OUTPUT "Frame\tX_ROI\tY_ROI\tAngle\tRevolution\r\n";
	## WHEN Refference file does exist
	if( $ref_exist eq 'yes' ) {
		print OUTPUT_RFD "Frame\tX_ROI\tY_ROI\r\n";
	}
	
	## wirte data
	$counter = 0;
	$skipped_lines = 0;
	while ( my $line = <INPUT> ) {
		$line =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
		($frame, $x_ROI, $y_ROI, $angle, $accum, $revolution) = ($1, $2, $3, $4, $5, $6);
		
		## write for non-drift corrected file.
		printf OUTPUT ("%5d\t%6.4f\t%6.4f\t%6.4f\t%6.4f\r\n",$frame,$x_ROI,$y_ROI,$angle,$revolution);
	
		## write for drift corrected file.
		if( $ref_exist eq 'yes' ) {
			$ref = <REF>;
			if(defined $ref) {
				$ref =~ /^([\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+).*/;
				$x_ROI = $x_ROI - $2 + 20; ## shift +20 for removing minus
				$y_ROI = $y_ROI - $3 + 20;
				printf OUTPUT_RFD ("%5d\t%6.4f\t%6.4f\r\n",$frame,$x_ROI,$y_ROI);
			}else{
				$skipped_lines++;
			}
		}
		$counter++;
	}
	if($skipped_lines) {
		$num_of_error++;
		push(@error_files, $file_Input);
		print"\tERROR: Refference is shorter than Input. Stoped at $counter. Skipped $skipped_lines lines.\n";
	}
	
	close(INPUT) || die "Cannot close: \n\t$!";
	close(OUTPUT) || die "Cannot close: \n\t$!";
	
	## for displaying result
	$file_Output =~ /([^\s^\/]+)$/;
	$file_name_only = $1;
	print "\tOUT >> \"$file_name_only\"\n";
	
	if( $ref_exist eq 'yes' ) {
		close(REF) || die "Cannot close: \n\t$!";
		close(OUTPUT_RFD) || die "Cannot close: \n\t$!";
		$num_of_rfd_file++;
		
		## for displaying result
		$file_Output_Rfd =~ /([^\s^\/]+)$/;
		$file_name_only = $1;
		print "\tOUT >> \"$file_name_only\"\n";
	}
}
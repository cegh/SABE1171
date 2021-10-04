#!/usr/bin/perl -w

# Version 1
# Author: Jaqueline Wang.
# Bioinformatician at Human Genome and Stem Cell Research Center - IB USP
# e-mail : jaqueytw@gmail.com
# Date : 12 November 2019

## CEGH USP VCF Filter
# Based on Tatiana O. M. S. and Michel N. script - CEGH_USP_VCF_Filter-v1.3.pl


### SCRIPT TO
# 1 - calculate and add CEGH filter
# 2 - correct the hemizygous variants
# 3 - deal with multiallelic sites
# 4 - use file after INCLUI COLUNAS, sample names are already included

use strict;

my ($file, $sex_file) = @ARGV;

if (!($file && $sex_file)) {
	die "Missing Input \
Usage : SABE_VCF_filter.pl <VCF_file> <SEX_file>\n"
}


### THE CODE

# Name the outfile
my $outfile = $file;
$outfile =~ s/\.txt$/.CEGHFilter.txt/;
open (OUT, ">$outfile") or die "Fail to open OUT file\n";


# Get the sample sex
my %sample_sex;
open (SEX, "$sex_file") or die "Fail to open SEX file\n";
my $head = <SEX>;

while (my $line = <SEX>) {
    chomp ($line);
    my @data = split (/\s/, $line);
    my $sample = $data[0];
    my $sex = $data[2];

    if (!exists($sample_sex{$sample})) {
	$sample_sex{$sample} = $sex;
    }

    else {
	print "Problem with $sample \n";
    }
}
close (SEX);


# Hg38 Pseudoautosomal regions in X and Y
my $pseudo_Y_10001 = 10000;
my $pseudo_Y_2649520 = 2781479;
my $pseudo_Y_59034050 = 56887902;
my $pseudo_Y_59363566 = 57217415;

my $pseudo_X_60001 = 10000;
my $pseudo_X_2699520 = 2781479;
my $pseudo_X_154931044 = 155701382;
my $pseudo_X_155260560 = 156030895;


# Declare variables:
# Countg means count genotypes
my ($count_samples, $count_variants, $number_columns);
my $countg_0_0 = 0; my $countg_0_1 = 0; my $countg_1_1 = 0;
my $countg_fd = 0; my $countg_fb = 0;
my $countg_fd_0_0 = 0; my $countg_fb_0_0 = 0;
my $countg_fd_0_1 = 0; my $countg_fb_0_1 = 0;
my $countg_fd_1_1 = 0; my $countg_fb_1_1 = 0;
my $countg_1_dot = 0; my $countg_dot_1 = 0; my $countg_0_dot  = 0; my $countg_dot_dot = 0; my $countg_dot_0 = 0; my $countg_other_options = 0;
my $countg_no_info = 0; # Count ./.:.:.:.
my $countg_homo = 0; # Count homo Alternative
my $countg_hemi = 0; # Count genotypes 1/1 in X-exclusive and Y-exclusive regions in Men
my $countg_hemi_0_0 = 0; # Count genotypes 0/0 in X-exclusive and Y-exclusive regions in Men
my $countg_homo_REF = 0; # Count homo Reference
my $countg_hetero = 0;
my $threshold_CEGH_filter = 1000; # MARILIA & MICHEL changed # Number of well genotyped alleles, choose based in your allele number
my $prop_bad_ALTg = 0; # Ratio of bad genotypes with alternate alleles (e.g. FDs- FBs-) over all genotypes with alternate alleles

# Coverage (depth) 
my $min_coverage_0_1 = 10;
my $min_coverage_homo = 7;
my $min_coverage_dot_1 = 5;

# Allele balance (ratios)
my $min_ac_ratio_0_1 = 0.3;
my $max_ac_ratio_0_1 = 0.7;
my $min_ac_ratio_0_0 = 0.9;
my $max_ac_ratio_0_0 = 1.1;
my $min_ac_ratio_1_1 = 0;
my $max_ac_ratio_1_1 = 0.1;

# Total coverage, frequencies, depth and allele count for ratio variables
my $sum; # Sum of reference and alternate reads inside a genotype call 
my $ac_ratio; # Ratio of reference reads over sum of reads inside a genotype call 
my $total_depth; # Total depth of every genotype with number of read info 
my $total_frequencies = 0; # Raw frequency before FD- and FB- flags, but after Male:Female info inclusion
my $total_frequencies_filtered = 0; # Clean frequency excluding FDs and FBs genotypes 
my $alleleALT_count_Raw = 0; # Number of alternate alleles in variant before filtering accounting for sex 
my $alleleTOT_number_Raw = 0; # Total number of alleles before Depth and Allele balance filtering and accounting for sex


# Open the in file and take de header to obtain some informations 
open (IN, "$file") or die "Fail to open IN file $file \n";
my $header = <IN>;
chomp ($header);
my @header_columns = split (/\t/, $header);


# Locate columns locations for parsing (ExonicFunc.refGene, GeneDetail.refGene, Func.refGene, AA.Change.refGene, FORMAT)
my $exon = 0; ++$exon until $header_columns[$exon] =~ /^ExonicFunc.refGene$/ or $exon > $#header_columns;++$exon;
my $gene_detail = 0; ++$gene_detail until $header_columns[$gene_detail] =~ /^GeneDetail.refGene$/ or $gene_detail > $#header_columns;++$gene_detail;
my $func = 0; ++$func until $header_columns[$func] =~ /^Func.refGene$/ or $func > $#header_columns;++$func;
my $AA_change = 0; ++$AA_change until $header_columns[$AA_change] =~ /^AAChange.refGene$/ or $AA_change > $#header_columns;++$AA_change;
my $format = 0; ++$format until $header_columns[$format] =~ /^FORMAT$/ or $format > $#header_columns;++$format;

# Relate the column number to the sample and the sex
my $count_col = 0;
my %sample_col_sex;
foreach my $key (@header_columns) {
	if (exists($sample_sex{$key})) {
		$sample_col_sex{$count_col} = $sample_sex{$key};
	}
	$count_col += 1;
}
$count_samples = $count_col - $format;

# Rename columns in header: {'Func.refGene'} for {'PredictedFunc.refGene'} and {'AA.Change.refGene'} for {'PredConsequence.refGene'}
my $f = $func - 1;
my $a = $AA_change - 1;
splice @header_columns, $f, 1, "PredictedFunc.refGene";
splice @header_columns, $a, 1, "PredConsequence.refGene";

# Remove {'ExonicFunc.refGene'} and {'GeneDetail.refGene'}
my $e = $exon - 1;
my $g = $gene_detail - 1;
splice @header_columns, $e, 1;
splice @header_columns, $g, 1;

# Insert new columns : Frequencies (raw and filtered), Allele count, Allele number, Heterozygous, Homozygous REF and ALT, Hemizygous counts and CEGH Filter before FORMAT column
my $frequencies_filtered_column = $format - 3;
my $frequencies_column = $format - 3;
my $alleleALT_count_Raw_column = $format - 3;
my $alleleTOT_number_Raw_column = $format - 3;
my $hetero_column = $format - 3;
my $homo_REF_column = $format - 3;
my $hemi_column = $format - 3;
my $homo_ALT_column = $format - 3;
my $CEGH_filter_column = $format - 3;
splice @header_columns, $frequencies_filtered_column, 0, "Frequencies_filtered";
splice @header_columns, $frequencies_column, 0, "Frequencies";
splice @header_columns, $alleleALT_count_Raw_column, 0, "Allele_ALT_count";
splice @header_columns, $alleleTOT_number_Raw_column, 0, "Allele_number";
splice @header_columns, $hetero_column, 0, "Heterozygous_count";
splice @header_columns, $homo_REF_column, 0, "HomozygousREF_count";
splice @header_columns, $hemi_column, 0, "Hemizygous_count";
splice @header_columns, $homo_ALT_column, 0, "HomozygousALT_count";
splice @header_columns, $CEGH_filter_column, 0, "CEGH_Filter";

$header = join ("\t", @header_columns);
print OUT "$header\n";


# Analyze each line at a time
while (my $line = <IN>) {
	chomp ($line);
	# Count number of variants
	++$count_variants;

	# Turn on Chr X or Y flag if in Chr X/Y variant, turn off if not
	my $flagX = "no";
	my $flagY = "no";

	if ($line =~ m/^ChrX/i or $line =~ m/^X/i) {
		$flagX = "yes";
	} 
	elsif ($line =~ m/^ChrY/i or $line =~ m/^Y/i){
		$flagY = "yes";
	}

	my @columns = split(/\t/,$line);

	# Grab start and end position of variant 
	my $start_position = $columns[1];
	my $end_position = $columns[2];

	# Declare variables for CEGH Filter flagging
	my $CEGH_filter_flag = "";

	# These are the counts of each genotype called in one given variant used for CEGH Filtering
	my $CEGH_filter_0_0 = 0;  my $CEGH_filter_0_1 = 0; my $CEGH_filter_1_1 = 0; 
	my $CEGH_filter_dot_1 = 0; my $CEGH_filter_1_dot = 0;
	my $CEGH_filter_fd_0_1 = 0; my $CEGH_filter_fb_0_1 = 0; 
	my $CEGH_filter_fd_0_0 = 0; my $CEGH_filter_fb_0_0 = 0;
	my $CEGH_filter_fd_dot_1 = 0;
	my $CEGH_filter_fb = 0; my $CEGH_filter_fd = 0; 
	my $CEGH_filter_well_genotyped = 0; # Count of genotypes with REFERENCE and ALTERNATE alleles that met filtering criteria, excluding X-exclusive and Y-exclusive regions (which count HEMI for men)
	my $CEGH_filter_well_genotyped_hemi_1_1 = 0; # Count of 1_1 genotypes in men in X- and Y-exclusive regions that met filtering criteria
	my $CEGH_filter_well_genotyped_hemi_0_0 = 0; # Count of 0_0 genotypes in men in X- and Y-exclusive regions that met filtering criteria
	my $CEGH_filter_fb_homo = 0; my $CEGH_filter_fd_homo = 0; 
	my $CEGH_filter_fb_hemi = 0; my $CEGH_filter_fd_hemi = 0;
	my $sum_1_1_0_1 = 0; # Sum of genotype counts containg alternate alleles 

	# specific genotypes counter - declaring all possibilities, including possibly spurious 0/1 in X- and Y-exclusive regions in Man
	my $n_0_0_good = 0; my $n_0_0_fd = 0; my $n_0_0_fb = 0; 
	my $n_0_1_good = 0; my $n_0_1_fd = 0; my $n_0_1_fb = 0; 
	my $n_1_1_good = 0; my $n_1_1_fd = 0; my $n_1_1_fb = 0;
	my $n_dot_1_good = 0; my $n_dot_1_fd = 0; 
	my $n_0_0_hemi_good = 0; my $n_0_0_hemi_fd = 0; my $n_0_0_hemi_fb = 0; 
	my $n_1_1_hemi_good = 0; my $n_1_1_hemi_fd = 0; my $n_1_1_hemi_fb = 0; 
	my $n_0_1_hemi_good; my $n_0_1_hemi_fb; my $n_0_1_hemi_fd;
	
	# Grab ExonicFunc.refGene column and Func.refGene column
	my %predicted_func; my $predicted_func_column;
	$predicted_func{"$columns[$f]"} = "$columns[$e]";

	# Grab GeneDetail.refGene column and AA.Change.refGene column
	my %pred_consequence; my $predicted_cons_column;
	$pred_consequence{"$columns[$a]"} = "$columns[$g]";

	my @exonic_func = values %predicted_func;
	my @func_ref = keys %predicted_func;		
	my @exonic_cons = values %pred_consequence;
	my @func_cons = keys %pred_consequence;

	# Check if key = exonic or not
	if ($func_ref[0] =~ m/^exonic$/) {
		$predicted_func_column = $exonic_func[0];
		$predicted_cons_column = $func_cons[0];
	} 
	
	else {
		$predicted_func_column = $func_ref[0];
		$predicted_cons_column = $exonic_cons[0];
	}


	# Add the flag for SEX and CEGH Filter (CF) on the format field
	$columns[$format - 1] = "$columns[$format - 1]:SEX:CF";

	#my @format_info = split (/:/, $columns[$format]);
        #my $GT = 0; ++$GT until $format_info[$GT] =~ /GT/ or $GT > $#format_info;
        #my $AD = 0; ++$AD until $format_info[$AD] =~ /AD/ or $AD > $#format_info;

	#For each line of the file, analyze one column at a time
	my $count_col = 0;
	foreach my $column (@columns) {
	
		# Grab alele count for reference and for alternative
		if ($column =~ m/^[01]\/[01]\:(\d+)\,(\d+)/) {
			my $ac_reference = $1;
			my $ac_alternative = $2;
			$sum = $ac_reference + $ac_alternative;
			$total_depth += $sum;	# Count for summary depth everytime there are reads

			# Attribute the sample sex
			my $sex = $sample_col_sex{$count_col};

			# Check coverage and ac_ratio for each genotype
			# First, hetero
			if ($column =~ m/^0\/1\:(\d+)\,(\d+)/) {
				++$countg_0_1; ++$CEGH_filter_0_1; 				
				
				# Count Homozygous and Hemizygous
				if ($flagX =~ m/no/i && $flagY =~ m/no/i) {
					++$countg_hetero;

					# If sum is lower than	min_coverage write FD
					if ($sum >= $min_coverage_0_1){
						$ac_ratio = $ac_reference / $sum;

						if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
							$columns[$count_col] = "$column:$sex:."; 
							++$CEGH_filter_well_genotyped; 
							++$sum_1_1_0_1; 
							++$n_0_1_good;
						} 

						else {
							$columns[$count_col] = "$column:$sex:FB"; 
							++$countg_fb_0_1; 
							++$CEGH_filter_fb; 
							++$CEGH_filter_fb_0_1; 
							++$n_0_1_fb;
						}
					} 

					else {
						$columns[$count_col] = "$column:$sex:FD"; 
						++$countg_fd_0_1; 
						++$CEGH_filter_fd; 
						++$CEGH_filter_fd_0_1; 
						++$n_0_1_fd;
					}
				} #if ($flagX =~ m/no/i && $flagY =~ m/no/i)

				elsif ($flagX =~ m/yes/i) {
					if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 
	
						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_0_1; 
							++$n_0_1_fd;
						}
					} #if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) 

					elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fb_0_1; ###AQUI
							++$n_0_1_fd;
						}
								
					} #elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

					elsif ($sex =~ /M/i) { # X Non-Pseudoautosomal region (X exclusive) FOR MAN
						++$countg_hemi;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1) {
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_hemi_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_hemi_fb;
							}
						} 
	
						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fb_0_1; ###AQUI
							++$n_0_1_hemi_fd;
						}
					} #elsif ($sex =~ /M/i)

					else { # X Non-Pseudoautosomal region (X exclusive) FOR WOMAN
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fb_0_1; ###AQUI
							++$n_0_1_fd;
						}
					} #else

				} #elsif ($flagX =~ m/yes/i)

				elsif ($flagY =~ m/yes/i) {

					if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_0_1; 
							++$n_0_1_fd;
						}
					} #if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)
				
					elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_0_1;
							++$n_0_1_fd;
						}
					} #elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)

					elsif ($sex =~ /M/i) { # Y Non-Pseudoautosomal region (Y exclusive) FOR MAN
						++$countg_hemi;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_hemi_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_hemi_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_0_1;
							++$n_0_1_hemi_fd;
						}
					} #elsif ($sex =~ /M/i)
					
					else { # Y Non-Pseudoautosomal region (Y exclusive) FOR WOMAN
						++$countg_hetero;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_0_1) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_1 && $ac_ratio <= $max_ac_ratio_0_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_0_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_0_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_0_1; 
								++$n_0_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_0_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_0_1;
							++$n_0_1_fd;
						}
					} #else

				} #elsif ($flagY =~ m/yes/i)

			} #if ($column =~ m/^0\/1\:(\d+)\,(\d+)/)

			elsif ($column =~ m/^1\/1\:(\d+)\,(\d+)/) {
				++$countg_1_1; ++$CEGH_filter_1_1;

				# Count Homozygous and Hemizygous
				if ($flagX =~ m/no/i && $flagY =~ m/no/i) {
					++$countg_homo;

					# If sum is lower than	min_coverage write FD
					if ($sum >= $min_coverage_homo) {
						$ac_ratio = $ac_reference / $sum;
						if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
							$columns[$count_col] = "$column:$sex:."; 
							++$CEGH_filter_well_genotyped; 
							++$sum_1_1_0_1; 
							++$n_1_1_good;
						} 

						else {
							$columns[$count_col] = "$column:$sex:FB"; 
							++$countg_fb_1_1; 
							++$CEGH_filter_fb; 
							++$CEGH_filter_fb_homo; 
							++$n_1_1_fb;
						}
					} 

					else {
						$columns[$count_col] = "$column:$sex:FD"; 
						++$countg_fd_1_1; 
						++$CEGH_filter_fd; 
						++$CEGH_filter_fd_homo; 
						++$n_1_1_fd;
					}
				} #if ($flagX =~ m/no/i && $flagY =~ m/no/i)

				elsif ($flagX =~ m/yes/i) {
			
					if ($sex =~ /M/i && $start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;				
						}
					} #if ($sex =~ /M/i && $start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520)

					elsif ($sex =~ /M/i && $start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1) {
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						}  
			
						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($sex =~ /M/i && $start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

					elsif ($sex =~ /M/i) { # X Non-Pseudoautosomal region (X exclusive)
						++$countg_hemi;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) { 
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped_hemi_1_1; 
								++$sum_1_1_0_1; 
								++$n_1_1_hemi_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_hemi; 
								++$n_1_1_hemi_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_hemi; 
							++$n_1_1_hemi_fd;
						}
					} #elsif ($sex =~ /M/i)

					elsif ($sex =~ /F/i && $start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col]= "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($sex =~ /F/i && $start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520)

					elsif ($sex =~ /F/i && $start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
						++$countg_homo;
						
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($sex =~ /F/i && $start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

					elsif ($sex =~ /F/i) { # X Non-Pseudoautosomal region (X exclusive)
						++$countg_homo;
						
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped;
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($sex =~ /F/i)

					elsif ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
						++$countg_homo; # No sex information

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520)

					elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
						++$countg_homo; # No sex information

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

					else { # X Non-Pseudoautosomal region (X exclusive) but no sex information
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) { 
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #else

				} #elsif ($flagX =~ m/yes/i)
				
				elsif ($flagY =~ m/yes/i) {

					if ($sex =~ /M/i && $start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #if ($sex =~ /M/i && $start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)
				
					elsif ($sex =~ /M/i && $start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
						++$countg_homo;
									
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1) {
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($sex =~ /M/i && $start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)

					elsif ($sex =~ /M/i) { # Y Non-Pseudoautosomal region (Y exclusive)
						++$countg_hemi;
												
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped_hemi_1_1; 
								++$sum_1_1_0_1; 
								++$n_1_1_hemi_good;
							} 

							else {
								$columns[$count_col]= "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_hemi; 
								++$n_1_1_hemi_fb;
							}

						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_hemi; 
							++$n_1_1_hemi_fd;
						}
					} #elsif ($sex =~ /M/i)

					elsif ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
						++$countg_homo;
							
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} #elsif ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)

					elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
						++$countg_homo;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1) { 
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$sum_1_1_0_1; 
								++$n_1_1_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_homo; 
								++$n_1_1_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_homo; 
							++$n_1_1_fd;
						}
					} # elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)

					else { # Y Non-Pseudoautosomal region (Y exclusive) but no sex information
						++$countg_hemi;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_1_1 && $ac_ratio <= $max_ac_ratio_1_1){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped_hemi_1_1; 
								++$sum_1_1_0_1; 
								++$n_1_1_hemi_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$countg_fb_1_1; 
								++$CEGH_filter_fb; 
								++$CEGH_filter_fb_hemi; 
								++$n_1_1_hemi_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$countg_fd_1_1; 
							++$CEGH_filter_fd; 
							++$CEGH_filter_fd_hemi; 
							++$n_1_1_hemi_fd;
						}
					} #else 

				} #elsif ($flagY =~ m/yes/i)
				
			} #elsif ($column =~ m/^1\/1\:(\d+)\,(\d+)/)
		
			# Last, homo 0/0			
			elsif ($column =~ m/^0\/0\:(\d+)\,(\d+)/) {
				++$countg_0_0; ++$CEGH_filter_0_0;

				if ($flagX =~ m/yes/i) {
					if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
						++$countg_homo_REF;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$n_0_0_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$CEGH_filter_fb_0_0; 
								++$n_0_0_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$CEGH_filter_fd_0_0; 
							++$n_0_0_fd;
						}									
					} #if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520)

					elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
						++$countg_homo_REF;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$n_0_0_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$CEGH_filter_fb_0_0; 
								++$n_0_0_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$CEGH_filter_fd_0_0; 
							++$n_0_0_fd;
						}									
					} #elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

					else {
						if ($sex =~ /M/i) {
							++$countg_hemi_0_0;

							# If sum is lower than	min_coverage write FD
							if ($sum >= $min_coverage_homo) {
								$ac_ratio = $ac_reference / $sum;
								if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
									$columns[$count_col] = "$column:$sex:."; 
									++$CEGH_filter_well_genotyped_hemi_0_0; 
									++$n_0_0_hemi_good;
								} 

								else {
									$columns[$count_col] = "$column:$sex:FB"; 
									++$CEGH_filter_fb_0_0; 
									++$n_0_0_hemi_fb;
								}
							} 

							else {
								$columns[$count_col] = "$column:$sex:FD"; 
								++$CEGH_filter_fd_0_0; 
								++$n_0_0_hemi_fd;
							}	
						} #if ($sex =~ /M/i) 

						else {
							++$countg_homo_REF;

							# If sum is lower than	min_coverage write FD
							if ($sum >= $min_coverage_homo) {
								$ac_ratio = $ac_reference / $sum;
								if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
									$columns[$count_col] = "$column:$sex:."; 
									++$CEGH_filter_well_genotyped; 
									++$n_0_0_good;
								} 

								else {
									$columns[$count_col] = "$column:$sex:FB"; 
									++$CEGH_filter_fb_0_0; 
									++$n_0_0_fb;
								}
							} 

							else {
								$columns[$count_col] = "$column:$sex:FD"; 
								++$CEGH_filter_fd_0_0; 
								++$n_0_0_fd;
							}	

						} #else 
					} #else
				
				} #if ($flagX =~ m/yes/i)

				elsif ($flagY =~ m/yes/i) {
					if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
						++$countg_homo_REF;

						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo) {
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$n_0_0_good;
							} 
							
							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$CEGH_filter_fb_0_0; 
								++$n_0_0_fb;
							}
						} 

						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$CEGH_filter_fd_0_0; 
							++$n_0_0_fd;
						}
					} #if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)

					elsif ($sex =~ /M/i && $start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
						++$countg_homo_REF;
	
						# If sum is lower than	min_coverage write FD
						if ($sum >= $min_coverage_homo){
							$ac_ratio = $ac_reference / $sum;
							if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
								$columns[$count_col] = "$column:$sex:."; 
								++$CEGH_filter_well_genotyped; 
								++$n_0_0_good;
							} 

							else {
								$columns[$count_col] = "$column:$sex:FB"; 
								++$CEGH_filter_fb_0_0; 
								++$n_0_0_fb;
							}
						} 
				
						else {
							$columns[$count_col] = "$column:$sex:FD"; 
							++$CEGH_filter_fd_0_0; 
							++$n_0_0_fd;
						}
					} #elsif ($sex =~ /M/i && $start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)
	
					else {
						if ($sex =~ /M/i) {
							++$countg_hemi_0_0;

							# If sum is lower than	min_coverage write FD
							if ($sum >= $min_coverage_homo){
								$ac_ratio = $ac_reference / $sum;
								if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
									$columns[$count_col] = "$column:$sex:."; 
									++$CEGH_filter_well_genotyped_hemi_0_0; 
									++$n_0_0_hemi_good;
								} 

								else {
									$columns[$count_col] = "$column:$sex:FB"; 
									++$CEGH_filter_fb_0_0; 
									++$n_0_0_hemi_fb;
								}
							} 

							else {
								$columns[$count_col] = "$column:$sex:FD"; 
								++$CEGH_filter_fd_0_0; 
								++$n_0_0_hemi_fd;
							}
						} #if ($sex =~ /M/i)

						else {
							++$countg_homo_REF;

							# If sum is lower than	min_coverage write FD
							if ($sum >= $min_coverage_homo){
								$ac_ratio = $ac_reference / $sum;
								if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
									$columns[$count_col] = "$column:$sex:."; 
									++$CEGH_filter_well_genotyped; 
									++$n_0_0_good;
								} 

								else {
									$columns[$count_col] = "$column:$sex:FB"; 
									++$CEGH_filter_fb_0_0; 
									++$n_0_0_fb;
								}
							} 

							else {
								$columns[$count_col] = "$column:$sex:FD"; 
								++$CEGH_filter_fd_0_0; 
								++$n_0_0_fd;
							}

						} #else
					} #else

				} #elsif ($flagY =~ m/yes/i)			

				else {
					++$countg_homo_REF;

					# If sum is lower than	min_coverage write FD
					if ($sum >= $min_coverage_homo){
						$ac_ratio = $ac_reference / $sum;
						if ($ac_ratio >= $min_ac_ratio_0_0 && $ac_ratio <= $max_ac_ratio_0_0){
							$columns[$count_col] = "$column:$sex:."; 
							++$CEGH_filter_well_genotyped; 
							++$n_0_0_good;
						} 

						else {
							$columns[$count_col] = "$column:$sex:FB"; 
							++$CEGH_filter_fb_0_0; 
							++$n_0_0_fb;
						}
					} 

					else {
						$columns[$count_col] = "$column:$sex:FD"; 
						++$CEGH_filter_fd_0_0; 
						++$n_0_0_fd;
					}
				} #else

			} #elsif ($column =~ m/^0\/0\:(\d+)\,(\d+)/)

		} #if ($column =~ m/^[01]\/[01]\:(\d+)\,(\d+)/)

		# If ./. or nothing
		else {

			# Attribute the sample sex
                        my $sex = $sample_col_sex{$count_col};
	
			# Is this line Multiallelic or poorly genotiped (no information)?
			if ($column =~ m/^[01\.]\/[01\.]:(\d+),(\d+)/) { #Multiallelic
				my $ac_reference = $1;
				my $ac_alternative = $2;
				$sum = $ac_reference + $ac_alternative;
				$total_depth += $sum;

				# Count Multiallelic options
				
				# Hetero with two alternate alleles
				if ($column =~ m/^1\/\.\:(\d+)\,(\d+)/) {
					++$countg_1_dot; ++$CEGH_filter_1_dot;

					# Count Homozygous and Hemizygous
					if ($flagX =~ m/no/i && $flagY =~ m/no/i) {
                                        	++$countg_hetero;

                                        	# If sum is lower than  min_coverage write FD
                                        	if ($sum >= $min_coverage_dot_1){
                                                	
							# Not possible to analyse the other allele
							# No FB
							$columns[$count_col] = "$column:$sex:MU";
                                                       	++$CEGH_filter_well_genotyped;
                                                       	++$sum_1_1_0_1;
							++$n_dot_1_good;
                                        	}

                                        	else {
                                                	$columns[$count_col] = "$column:$sex:MU-FD";
							#++$countg_fd_0_1;
                                                	++$CEGH_filter_fd;
                                                	++$CEGH_filter_fd_dot_1;
							++$n_dot_1_fd;
                                        	}
                                	} #if ($flagX =~ m/no/i && $flagY =~ m/no/i)
 
					elsif ($flagX =~ m/yes/i) {

						if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
        	                                        ++$countg_hetero;
	
                	                                # If sum is lower than  min_coverage write FD
                        	                        if ($sum >= $min_coverage_dot_1){

								# Not possible to analyse the other allele
								# No FB
                                                	        $columns[$count_col] = "$column:$sex:MU"; 
                                                        	++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
								++$n_dot_1_good;
                                        	        } 
        
                                                	else {
                                                        	$columns[$count_col] = "$column:$sex:MU-FD"; 
								#++$countg_fd_0_1; 
        	                                                ++$CEGH_filter_fd; 
                	                                        ++$CEGH_filter_fd_dot_1; 
								++$n_dot_1_fd;
                                	                }
                                        	} #if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) 

						elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
	                                                ++$countg_hetero;
	
        	                                        # If sum is lower than  min_coverage write FD
                	                                if ($sum >= $min_coverage_dot_1){
                                        	
								# Not possible to analyse the other allele
								# No FB
								$columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                       	        ++$sum_1_1_0_1; 
                                                               	++$n_dot_1_good;
                                                	} 
	
        	                                        else {
                	                                        $columns[$count_col] = "$column:$sex:MU-FD"; 
                        	                                #++$countg_fd_0_1; 
                                	                        ++$CEGH_filter_fd; 
                                        	                ++$CEGH_filter_fd_dot_1; 
                                                	        ++$n_dot_1_fd;
	                                                }
        	                                                        
                	                        } #elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

						elsif ($sex =~ /M/i) { # X Non-Pseudoautosomal region (X exclusive) FOR MAN
        	                                        ++$countg_hemi;

                	                                # If sum is lower than  min_coverage write FD
                        	                        if ($sum >= $min_coverage_dot_1) {

								# Not possible to analyse the other allele
								# No FB
	                                               		$columns[$count_col] = "$column:$sex:MU"; 
                                                       	        ++$CEGH_filter_well_genotyped; 
                                                               	++$sum_1_1_0_1; 
                                                                #++$n_0_1_hemi_good;
        	                                        } 
        
                	                                else {
                        	                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                	                        #++$countg_fd_0_1; 
                                        	                ++$CEGH_filter_fd; 
                                                	        ++$CEGH_filter_fd_dot_1; 
                                                        	#++$n_0_1_hemi_fd;
	                                                }
        	                                } #elsif ($sex =~ /M/i)

						else { # X Non-Pseudoautosomal region (X exclusive) FOR WOMAN
	                                                ++$countg_hetero;

        	                                        # If sum is lower than  min_coverage write FD
                	                                if ($sum >= $min_coverage_dot_1) {

								# Not possible to analyse the other allele
								# No FB
                                       	                        $columns[$count_col] = "$column:$sex:MU"; 
                                               	                ++$CEGH_filter_well_genotyped; 
                                                       	        ++$sum_1_1_0_1; 
                                                               	++$n_dot_1_good;
	                                                } 

        	                                        else {
                	                                        $columns[$count_col] = "$column:$sex:MU-FD"; 
                        	                                #++$countg_fd_0_1; 
                                	                        ++$CEGH_filter_fd; 
                                        	                ++$CEGH_filter_fd_dot_1; 
                                                	        ++$n_dot_1_fd;
	                                                }
        	                                } #else

					} #elsif ($flagX =~ m/yes/i)

					elsif ($flagY =~ m/yes/i) {
					
						if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
	                                                ++$countg_hetero;

        	                                        # If sum is lower than  min_coverage write FD
                	                                if ($sum >= $min_coverage_dot_1) {
	
								# Not possible to analyse the other allele
								# No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                	} 

	                                                else {
        	                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                	                                        #++$countg_fd_0_1; 
                        	                                ++$CEGH_filter_fd; 
                                	                        ++$CEGH_filter_fd_dot_1; 
                                        	                ++$n_dot_1_fd;
                                                	}
	                                        } #if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)

						elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
	                                                ++$countg_hetero;

        	                                        # If sum is lower than  min_coverage write FD
                	                                if ($sum >= $min_coverage_dot_1) {
								# Not possible to analyse the other allele
								# No FB							
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                	} 

	                                                else {
        	                                                $columns[$count_col] = "$column:$sex:FD-MU"; 
                	                                        #++$countg_fd_0_1; 
                        	                                ++$CEGH_filter_fd; 
                                	                        ++$CEGH_filter_fd_0_1; 
                                        	                ++$n_dot_1_fd;
                                                	}
	                                        } #elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)

						elsif ($sex =~ /M/i) { # Y Non-Pseudoautosomal region (Y exclusive) FOR MAN
	                                                ++$countg_hemi;

        	                                        # If sum is lower than  min_coverage write FD
        	                                        if ($sum >= $min_coverage_dot_1){
								# Not possible to analyse the other allele
								# No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                #++$n_0_1_hemi_good;
                                                        }

	                                                else {
        	                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                	                                        #++$countg_fd_0_1; 
                        	                                ++$CEGH_filter_fd; 
                                	                        ++$CEGH_filter_fd_0_1; 
                                        	                #++$n_0_1_hemi_fd;
                                                	}
	                                        } #elsif ($sex =~ /M/i)

						else { # Y Non-Pseudoautosomal region (Y exclusive) FOR WOMAN
	                                                ++$countg_hetero;

        	                                        # If sum is lower than  min_coverage write FD
                	                                if ($sum >= $min_coverage_0_1) {
								# Not possible to analyse the other allele
								# No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        }

	                                                else {
        	                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                	                                        #++$countg_fd_0_1; 
                        	                                ++$CEGH_filter_fd; 
                                	                        ++$CEGH_filter_fd_0_1; 
                                        	                ++$n_dot_1_fd;
                                                	}
	                                        } #else

					} # elsif ($flagY =~ m/yes/i)

				} #elsif ($column =~ m/^1\/\.\:(\d+)\,(\d+)/)

				elsif ($column =~ m/^\.\/1\:(\d+)\,(\d+)/) {
					++$countg_dot_1; ++$CEGH_filter_dot_1;

					# Count Homozygous and Hemizygous
                                        if ($flagX =~ m/no/i && $flagY =~ m/no/i) {
                                                ++$countg_hetero;

                                                # If sum is lower than  min_coverage write FD
                                                if ($sum >= $min_coverage_dot_1){
                                                        
                                                        # Not possible to analyse the other allele
                                                        # No FB
                                                        $columns[$count_col] = "$column:$sex:MU";
                                                        ++$CEGH_filter_well_genotyped;
                                                        ++$sum_1_1_0_1;
                                                        ++$n_dot_1_good;
                                                }

                                                else {
                                                        $columns[$count_col] = "$column:$sex:MU-FD";
                                                        #++$countg_fd_0_1;
                                                        ++$CEGH_filter_fd;
                                                        ++$CEGH_filter_fd_dot_1;
                                                        ++$n_dot_1_fd;
                                                }
                                        } #if ($flagX =~ m/no/i && $flagY =~ m/no/i)

					elsif ($flagX =~ m/yes/i) {
				
						if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) { # X Pseudoautosomal region (X PAR1)
                                                        ++$countg_hetero;
        
                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1){

                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        } 
        
                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_dot_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                } #if ($start_position >= $pseudo_X_60001 && $end_position <= $pseudo_X_2699520) 

						elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560) { # X Pseudoautosomal region (X PAR2)
                                                        ++$countg_hetero;
        
                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1){
                                                
                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        } 
        
                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_dot_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                                        
                                                } #elsif ($start_position >= $pseudo_X_154931044 && $end_position <= $pseudo_X_155260560)

						elsif ($sex =~ /M/i) { # X Non-Pseudoautosomal region (X exclusive) FOR MAN
                                                        ++$countg_hemi;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1) {

                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                #++$n_0_1_hemi_good;
                                                        } 
        
                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_dot_1; 
                                                                #++$n_0_1_hemi_fd;
                                                        }
                                                } #elsif ($sex =~ /M/i)

						else { # X Non-Pseudoautosomal region (X exclusive) FOR WOMAN
                                                        ++$countg_hetero;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1) {

                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        } 

                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_dot_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                } #else


					} #elsif ($flagX =~ m/yes/i)

					elsif ($flagY =~ m/yes/i) {
	
						if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520) { # Y Pseudoautosomal region (Y PAR1)
                                                        ++$countg_hetero;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1) {
        
                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        } 

                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_dot_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                } #if ($start_position >= $pseudo_Y_10001 && $end_position <= $pseudo_Y_2649520)

						elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566) { # Y Pseudoautosomal region (Y PAR2)
                                                        ++$countg_hetero;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1) {
                                                                # Not possible to analyse the other allele
                                                                # No FB                                                 
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        } 

                                                        else {
                                                                $columns[$count_col] = "$column:$sex:FD-MU"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_0_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                } #elsif ($start_position >= $pseudo_Y_59034050 && $end_position <= $pseudo_Y_59363566)

						elsif ($sex =~ /M/i) { # Y Non-Pseudoautosomal region (Y exclusive) FOR MAN
                                                        ++$countg_hemi;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_dot_1){
                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                #++$n_0_1_hemi_good;
                                                        }

                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_0_1; 
                                                                #++$n_0_1_hemi_fd;
                                                        }
                                                } #elsif ($sex =~ /M/i)

						else { # Y Non-Pseudoautosomal region (Y exclusive) FOR WOMAN
                                                        ++$countg_hetero;

                                                        # If sum is lower than  min_coverage write FD
                                                        if ($sum >= $min_coverage_0_1) {
                                                                # Not possible to analyse the other allele
                                                                # No FB
                                                                $columns[$count_col] = "$column:$sex:MU"; 
                                                                ++$CEGH_filter_well_genotyped; 
                                                                ++$sum_1_1_0_1; 
                                                                ++$n_dot_1_good;
                                                        }

                                                        else {
                                                                $columns[$count_col] = "$column:$sex:MU-FD"; 
                                                                #++$countg_fd_0_1; 
                                                                ++$CEGH_filter_fd; 
                                                                ++$CEGH_filter_fd_0_1; 
                                                                ++$n_dot_1_fd;
                                                        }
                                                } #else

					} #elsif ($flagY =~ m/yes/i)

				} #elsif ($column =~ m/^\.\/1\:(\d+)\,(\d+)/)

				# MU-REP - Genotypes already analysed in another position
				elsif ($column =~ m/^[0]\/[\.]\:(\d+)\,(\d+)/){
					++$countg_0_dot;
					$columns[$count_col] = "$column:$sex:MU-REP";
				} 
				elsif ($column =~ m/^[\.]\/[0]\:(\d+)\,(\d+)/){
					++$countg_dot_0;
					$columns[$count_col] = "$column:$sex:MU-REP";
				} 
				elsif  ($column =~ m/^[\.]\/[\.]\:(\d+)\,(\d+)/){
					++$countg_dot_dot;
					$columns[$count_col] = "$column:$sex:MU-REP";
				}
				else {
					++$countg_other_options;
					$columns[$count_col] = "$column:$sex:MU-OTH";
				}
			}

			else {
				# Count Poorly genotyped 
				if ($column =~ m/^[\.]\/[\.]\:(\.)\:(\.)/) {
					$columns[$count_col] = "$column:$sex:PO";
					++$countg_no_info;
				}
			}
		} #else

		$count_col += 1;
	}

	# CEGH Filter flagging
	my $sum_total = $CEGH_filter_1_1 + $CEGH_filter_0_1 + $CEGH_filter_1_dot + $CEGH_filter_dot_1; # Only count 0/1 and 1/1 genotypes before flagging
	my $sum_FB_FD =  $CEGH_filter_fd + $CEGH_filter_fb; # Only count 0/1 and 1/1 genotypes flagged with FD- and FB-
	my $CEGH_filter_well_genotyped_alleles = $CEGH_filter_well_genotyped + $CEGH_filter_well_genotyped + $CEGH_filter_well_genotyped_hemi_0_0 + $CEGH_filter_well_genotyped_hemi_1_1;
	my $sum_good_and_bad_altg = $sum_1_1_0_1 + $CEGH_filter_well_genotyped_hemi_1_1 + $sum_FB_FD; 

	# Forcing denominator to be different from 0
	if ( $sum_good_and_bad_altg == 0) {
		$prop_bad_ALTg = 0;
	} 

	else {
		$prop_bad_ALTg = $sum_FB_FD / $sum_good_and_bad_altg;
	}

	if ( $prop_bad_ALTg == 1) {
		if ($CEGH_filter_fd > 0.5 * $sum_FB_FD) { 
			$CEGH_filter_flag = "FDP";
		} 

		else { # If in the rare occasion 50% of FD and FB happens, flag it with FAB 
			$CEGH_filter_flag = "FAB";
		}

	} 

	elsif ($prop_bad_ALTg >= 0) {
		if ($sum_1_1_0_1 > 0 && $sum_FB_FD > 0.5 * $sum_total && $CEGH_filter_well_genotyped_alleles >= $threshold_CEGH_filter) {
			$CEGH_filter_flag = "WK";
		} 

		elsif ($sum_1_1_0_1 > 0 && $sum_FB_FD > 0.5 * $sum_total && $CEGH_filter_well_genotyped_alleles < $threshold_CEGH_filter) {
			$CEGH_filter_flag = "WK-LowCall";
		} 

		elsif ($sum_1_1_0_1 > 0 && $sum_FB_FD > 0.1 * $sum_total && $sum_FB_FD <= 0.5 * $sum_total && $CEGH_filter_well_genotyped_alleles >= $threshold_CEGH_filter) {
			$CEGH_filter_flag = "SR";
		} 

		elsif ($sum_1_1_0_1 > 0 && $sum_FB_FD <= 0.1 * $sum_total && $CEGH_filter_well_genotyped_alleles >= $threshold_CEGH_filter) {
			$CEGH_filter_flag = "vSR";
		} 

		elsif ($sum_1_1_0_1 > 0 && $sum_FB_FD > 0.1 * $sum_total && $sum_FB_FD <= 0.5 * $sum_total && $CEGH_filter_well_genotyped_alleles < $threshold_CEGH_filter) {
			$CEGH_filter_flag = "SR-LowCall";
		} 

		elsif ($sum_1_1_0_1 > 0 && $sum_FB_FD <= 0.1 * $sum_total && $CEGH_filter_well_genotyped_alleles < $threshold_CEGH_filter) {
			$CEGH_filter_flag = "vSR-LowCall";
		}
	}

	# Calculate total frequencies and allele count/number
	my $sum_hetero_filtered = $countg_hetero - $CEGH_filter_fd_0_1 - $CEGH_filter_fb_0_1; # Number of hetero not flagged with FD- and FB-
	my $sum_hemi_filtered = $countg_hemi - $CEGH_filter_fd_hemi - $CEGH_filter_fb_hemi; # Number of hemi not flagged with FD- and FB-
	my $sum_homo_filtered = $countg_homo - $CEGH_filter_fd_homo - $CEGH_filter_fb_homo; # Number of homo alt not flagged with FD- and FB-
	my $sum_homo_REF_flagged = $CEGH_filter_fb_0_0 + $CEGH_filter_fd_0_0; # Number of homo ref flagged with FD- and FB-
	my $sum_homo_REF_filtered = $countg_homo_REF - $sum_homo_REF_flagged; # Number of homo ref not flagged with FD- and FB-
	my $sum_up = $countg_hetero + $countg_hemi + $countg_homo + $countg_homo;
	my $sum_up_filtered = $sum_hetero_filtered + $sum_hemi_filtered + $sum_homo_filtered + $sum_homo_filtered;
	my $sum_down = $countg_homo + $countg_homo + $countg_homo_REF + $countg_homo_REF + $countg_hetero + $countg_hetero + $countg_hemi;
	my $sum_down_filtered = $sum_homo_filtered + $sum_homo_filtered + $sum_homo_REF_filtered + $sum_homo_REF_filtered + $sum_hetero_filtered + $sum_hetero_filtered + $sum_hemi_filtered;		

	#frequencies with specific genotypes counter
	#experimental allele number raw
	my $allraw = $n_0_0_good + $n_0_0_fd + $n_0_0_fb + $n_0_0_good + $n_0_0_fd + $n_0_0_fb + $n_0_1_good + $n_0_1_fd + $n_0_1_fb + $n_0_1_good + $n_0_1_fd + $n_0_1_fb + $n_1_1_good + $n_1_1_fd + $n_1_1_fb + $n_1_1_good + $n_1_1_fd + $n_1_1_fb + $n_0_0_hemi_good + $n_0_0_hemi_fd + $n_0_0_hemi_fb + $n_1_1_hemi_good + $n_1_1_hemi_fd + $n_1_1_hemi_fb;
	my $alternatesraw = $n_0_1_good + $n_0_1_fd + $n_0_1_fb + $n_1_1_good + $n_1_1_fd + $n_1_1_fb + $n_1_1_good + $n_1_1_fd + $n_1_1_fb + $n_1_1_hemi_good + $n_1_1_hemi_fd + $n_1_1_hemi_fb;

	if ($allraw == 0) { 
		$total_frequencies = 0;
	} 

	else {
	    	$total_frequencies = $alternatesraw / $allraw;
	}

	my $allgood = 2*($n_0_0_good) + 2*($n_0_1_good) + 2*($n_1_1_good) + 1*($n_0_0_hemi_good) + 1*($n_1_1_hemi_good);
	my $alternatesgood = $n_0_1_good + 2*($n_1_1_good) + $n_1_1_hemi_good;
	if ($allgood == 0) { 
		$total_frequencies_filtered = 0;
	} 
	
	else {
		$total_frequencies_filtered = $alternatesgood / $allgood;
	}

	# Atualizando as contagens de hemizigotos também
	my $hemi_raw_counts = $n_1_1_hemi_good + $n_1_1_hemi_fd + $n_1_1_hemi_fb;
	
#=pod

#	if ($sum_down == 0) { 
#		$total_frequencies = 0;
#	} 

#	else {
#	    	$total_frequencies = $sum_up / $sum_down;
#	}

#	if ($sum_down_filtered == 0) { 
#		$total_frequencies_filtered = 0;
#	} 
	
#	else {
#		$total_frequencies_filtered = $sum_up_filtered / $sum_down_filtered;
#	}

#=cut

	$alleleALT_count_Raw = $countg_hetero + $countg_hemi + $countg_homo + $countg_homo;
	$alleleTOT_number_Raw = $countg_hetero + $countg_hetero + $countg_hemi + $countg_homo + $countg_homo + $countg_homo_REF + $countg_homo_REF;

	# Insert PredictedFunc.refGene and PredConsequence.refGene new columns
	splice @columns, $f, 1, $predicted_func_column;
	splice @columns, $a, 1, $predicted_cons_column;

	# Remove ExonicFunc.refGene and GeneDetail.refGene columns
	splice @columns, $e, 1;
	splice @columns, $g, 1;

        # Round frequencies RAW and filtered to 6 digits after decimal point
        my $total_frequencies_filtered_6 = sprintf("%.6f", $total_frequencies_filtered);
        my $total_frequencies_6 = sprintf("%.6f", $total_frequencies);


	# Insert Frequencies, allele count/number, Hetero/Homo/Hemizygous count column and Filter-in-house column in the line of our filtered outfile
	splice @columns, $frequencies_filtered_column, 0, "$total_frequencies_filtered_6"; $total_frequencies_filtered = 0;
	splice @columns, $frequencies_column, 0, "$total_frequencies_6"; $total_frequencies = 0;
	splice @columns, $alleleALT_count_Raw_column, 0, "$alternatesraw"; $alternatesraw = 0;
	splice @columns, $alleleTOT_number_Raw_column, 0, "$allraw"; $allraw = 0; #experimental aqui
	splice @columns, $hetero_column, 0, "$countg_hetero"; $countg_hetero = 0;
	splice @columns, $homo_REF_column, 0, "$countg_homo_REF"; $countg_homo_REF = 0; # Count homo reference (0/0)
	splice @columns, $hemi_column, 0, "$hemi_raw_counts"; $hemi_raw_counts = 0;
	splice @columns, $homo_ALT_column, 0, "$countg_homo"; $countg_homo = 0; # Count homo alternative (1/1)
	splice @columns, $CEGH_filter_column, 0, "$CEGH_filter_flag";
	my $filtered_line = join ("\t", @columns);

	print OUT "$filtered_line\n";
}

close (IN);

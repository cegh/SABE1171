#!/usr/bin/perl -w

use strict;

my ($lof, $vqsr, $vcf) = @ARGV;

if (!($lof && $vqsr && $vcf)) {
	die "Missing Input \
Usage: Pre_CEGHFilter_Recalculation_Loftee_VQSR-AS.pl <LOFTEE_file> <VQSR-AS_file> <IN_VCF> \n"
}

### THE CODE

### Open LOFTEE filtered file and keep annotations
open (LOF, "$lof") or die "Fail to open $lof file\n";
my %loft;
<LOF>;

while (my $line = <LOF>) {
	chomp ($line);
	my @data = split (/\t/, $line);
	
	if (!exists($loft{$data[0]})) {
		$loft{$data[0]}{ feat } = $data[1];
		$loft{$data[0]}{ cons } = $data[2];
		$loft{$data[0]}{ cdna } = $data[3];
		$loft{$data[0]}{ lof1 } = $data[4];
		$loft{$data[0]}{ lof2 } = $data[5];
		$loft{$data[0]}{ lof3 } = $data[6];
		$loft{$data[0]}{ lof4 } = $data[7];
	}
	else {
		print "Repetiu $data[0] em LOFTEE\n";
	}
}
close (LOF);


### Open VQSR-AS filtered file and keep information
open (VQSR, "$vqsr") or die "Fail to open $vqsr file \n";
my %vqsr;
<VQSR>;

while (my $line = <VQSR>) {
	chomp ($line);
	my @data = split (/\t/, $line);
	my $cppra = "$data[0]:$data[1]-$data[2]:$data[3]/$data[4]";

	if (!exists($vqsr{$cppra})) {
		$vqsr{$cppra} = $data[5];
	}	
	else {
		$vqsr{$cppra} = "$vqsr{$cppra};$data[5]";
	}
}
close (VQSR);


# Open VCF file
open (VCF, "$vcf") or die "Fail to open $vcf file \n";

my $out1 = $vcf;
$out1 =~ s/\.txt$/.Recalculation_Loftee_VQSR-AS.txt/;
open (OUT1, ">$out1") or die "Fail to open out file\n";

my $out2 = $vcf;
$out2 =~ s/\.txt$/.Sample187547917.txt/;
open (OUT2, ">$out2") or die "Fail to open sample 187547917 file\n";

my $chr = 0; my $ini = 0; my $end = 0; my $ref = 0; my $alt = 0;
my $cac = 0; my $caf = 0; my $can = 0; my $fil = 0; my $pli = 0;
my $sam_cont = 0;
my $cont_lof = 0;

while (my $line = <VCF>) {
	chomp ($line);
	my @data = split (/\t/, $line);

	if ($line =~ m/^Chr/) {
		++$chr until $data[$chr] =~ m/^Chr$/ or $chr > $#data;
		++$ini until $data[$ini] =~ m/^Start$/ or $ini > $#data;
		++$end until $data[$end] =~ m/^End$/ or $end > $#data;
		++$ref until $data[$ref] =~ m/^Ref$/ or $ref > $#data;
		++$alt until $data[$alt] =~ m/^Alt$/ or $alt > $#data;
		++$fil until $data[$fil] =~ m/^FILTER$/ or $fil > $#data;
		++$sam_cont until $data[$sam_cont] =~ m/^187547917$/ or $sam_cont > $#data;
		++$cac until $data[$cac] =~ m/^AC$/ or $cac > $#data;
		++$caf until $data[$caf] =~ m/^AF$/ or $caf > $#data;
		++$can until $data[$can] =~ m/^AN$/ or $can > $#data;
		++$pli until $data[$pli] =~ m/^pLI$/ or $pli > $#data;

		splice (@data, $sam_cont, 1);
		splice (@data, $cac, 0, "FILTER_VQSR-AS");
		splice (@data, $pli, 0, "Loftee_feature", "Loftee_Consequence", "Loftee_cDNA_position", "Loftee_LoF", "Loftee_LoF_filter", "Loftee_LoF_flags", "Loftee_LoF_info");

		my $new = join ("\t", @data);
		print OUT1 "$new\n"; 
		print OUT2 "$line\n";
	}
	else {

		my $cppra = "$data[$chr]:$data[$ini]-$data[$end]:$data[$ref]/$data[$alt]";
		my @newdata = @data;

		my @sam = split (/:/, $newdata[$sam_cont]);
		if ($sam[0] =~ m/\/1$/) {
			if ($sam[0] eq "0/1") {
				$newdata[$cac] = $newdata[$cac] - 1;
				$newdata[$can] = $newdata[$can] - 2;
				if ($newdata[$cac] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];

					splice (@newdata, $sam_cont, 1);
				}
				else {
					print OUT2 "$line\n";
					next;
				}				
			}
			elsif ($sam[0] eq "1/1") {
				$newdata[$cac] = $newdata[$cac] - 2;
				$newdata[$can] = $newdata[$can] - 2;
				if ($newdata[$cac] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];

					splice (@newdata, $sam_cont, 1);
				}
				else {
					print OUT2 "$line\n";
					next;
				}
			}
			elsif ($sam[0] eq "./1") {
				$newdata[$cac] = $newdata[$cac] - 1;
				$newdata[$can] = $newdata[$can] - 1;
				if ($newdata[$cac] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];

					splice (@newdata, $sam_cont, 1);
				}
				else {
					print OUT2 "$line\n";
					next;
				}
			}
			else {
				print "$cppra\t$data[$sam_cont]\n";
				next;
			}
		}
		else {
			if ($sam[0] eq "0/0") {
				$newdata[$can] = $newdata[$can] - 2;
				if ($newdata[$can] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];

					splice (@newdata, $sam_cont, 1);
				}
				else {
					print OUT2 "$line\n";
					next;
				}
			}
			elsif (($sam[0] eq "./0") || ($sam[0] eq "0/.")) {
				$newdata[$can] = $newdata[$can] - 1;
				if ($newdata[$can] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];

					splice (@newdata, $sam_cont, 1);
				}
				else {
					print OUT2 "$line\n";
					next;
				}
			}
			elsif ($sam[0] eq "./.") {
				if ($newdata[$cac] > 0) {
					$newdata[$caf] = $newdata[$cac]/$newdata[$can];
					$newdata[$caf] = sprintf '%.6f', $newdata[$caf];
	
					splice (@newdata, $sam_cont, 1);

				}
				else {
					print OUT2 "$line\n";
					next;
				}
			}
			else {
				print "$cppra\t$data[$sam_cont]\n";
				next;
			}
		}

		splice (@newdata, $cac, 0, $vqsr{$cppra});

		if (exists($loft{$cppra})) {
			splice (@newdata, $pli, 0, $loft{$cppra}{ feat }, $loft{$cppra}{ cons }, $loft{$cppra}{ cdna }, $loft{$cppra}{ lof1 }, $loft{$cppra}{ lof2 }, $loft{$cppra}{ lof3 }, $loft{$cppra}{ lof4 });
			$cont_lof += 1;
		}
		else {
			splice (@newdata, $pli, 0, "", "", "", "", "", "", "");
		}

		my $new = join ("\t", @newdata);
		print OUT1 "$new\n";

	}	
}
close (VCF);
close (OUT1);
close (OUT2);

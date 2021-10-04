#!/usr/bin/perl

# Version 1
# Author : Jaqueline Wang
# Bioinformatician at Human Genome and Stem Cell Research Center - IB USP
# e-mail : jaqueytw@gmail.com
# Date : 15 September 2020

### SCRIPT TO: Split multiallelic rows, including INFO column for dbSNP

# Observation
# 1 - INFO: AC, AN and AF and others
# 2 - FORMAT: GT, AD, DP, PL, PGT and PID

use strict;

my ($arquivo) = @ARGV; 

if (!($arquivo)) {
    die "Missing input \
Usage : Split_multiallelic_Haplotype_SABE-dbSNP.pl <VCF_file> \n";
}

### THE CODE

my $out_file = $arquivo;
$out_file =~ s/.vcf/.split.vcf/;

open (IN, "$arquivo") or die "Fail to open IN file $arquivo \n";
open (OUT,">$out_file") or die "Fail to open OUT file $out_file \n";

my $totalLinePre = 0;
my $totalLinePos = 0;
my $totalLineSplit = 0;
my $totalLineSingle = 0;


my $alt = 0;
my $inf = 0;
my $for = 0;

while (my $line = <IN>) {

	chomp ($line);

	### Keep columns info
	if ($line =~ m/#CHROM/) {
		print OUT "$line\n";
		my @dados = split(/\t/, $line);
		
		++$for until $dados[$for] =~ m/^FORMAT$/ or $for > $#dados;
		++$inf until $dados[$inf] =~ m/^INFO$/ or $inf > $#dados;
		++$alt until $dados[$alt] =~ m/^ALT$/ or $alt > $#dados;
	}

	### Print other lines of VCF header
	elsif ($line =~ m/##/) {
		print OUT "$line\n";
	} #elsif ($line =~ m/##/) 


	### Variants 
	elsif ($line =~ m/^chr/) {
        	$totalLinePre += 1; 

		my @dados = split(/\t/, $line);
		my @alt = split(/,/, $dados[$alt]);
		my $num_alt = scalar(@alt);

		my @format = split(/:/, $dados[$for]);
		my $count_form = 0;
		my $exist_PGT = 0;

		my ($AD_col, $DP_col, $GQ_col, $PL_col, $SAC_col, $PGT_col, $PID_col);

		### For each line, keep the FORMAT annotations
		foreach my $form_info (@format) {

			if ($form_info eq "AD") {
				$AD_col = $count_form;
	    		}
	    		elsif ($form_info eq "DP") {
			    	$DP_col = $count_form;
	    		}
			elsif ($form_info eq "GQ") {
				$GQ_col = $count_form;
	    		}
	    		elsif ($form_info eq "PL") {
				$PL_col = $count_form;
	    		}
	    		elsif ($form_info eq "SAC") {        
                		$SAC_col = $count_form;
            		}
	    		elsif ($form_info eq "PGT") {
				$PGT_col = $count_form;
				$exist_PGT += 1;
	    		}
	    		elsif ($form_info eq "PID") {
				$PID_col = $count_form;
	    		}
		
	    		$count_form += 1;
		}


		### IF PGT does NOT exists
		if ($exist_PGT == 0) {
	    		$dados[$for] = "GT:AD:DP:GQ:PL:SAC";

			### IF there is more than ONE ALT
	    		if ($num_alt > 1) {
				$totalLineSplit += 1;



				my @info = split(/;/, $dados[$inf]);

				my (@AC, @AF);
				my (@AS_Bas, @AS_FS, @AS_Fil, @AS_Inb, @AS_MQ, @AS_MQR, @AS_QD, @AS_Rea, @AS_SOR, @AS_VQS, @AS_cul);
				my (@MLEAC, @MLEAF);
				my ($ac_col, $af_col);
				my ($as_bas, $as_fs, $as_fil, $as_inb, $as_mq, $as_mqr, $as_qd, $as_rea, $as_sor, $as_vqs, $as_cul);
				my ($mleac_col, $mleaf_col);

				my $count_info = 0;
				foreach my $info_data (@info) {

		    			if ($info_data =~ m/^AC=/) {
						my $ac = $info_data;
						$ac =~ s/^AC=//;
						@AC = split (/,/, $ac);
						$ac_col = $count_info;
		    			}
		    			elsif ($info_data =~ m/^AF=/) {
						my $af = $info_data;
						$af =~ s/^AF=//;
						@AF = split (/,/, $af);
						$af_col = $count_info;
		    			}
					elsif ($info_data =~ m/^AS_BaseQRankSum=/) {
						my $as = $info_data;
						$as =~ s/^AS_BaseQRankSum=//;
						@AS_Bas = split (/,/, $as);
						$as_bas = $count_info;
					}		
					elsif ($info_data =~ m/^AS_FS=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_FS=//;
                                                @AS_FS = split (/,/, $as);
                                                $as_fs = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_FilterStatus=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_FilterStatus=//;
                                                @AS_Fil = split (/,/, $as);
                                                $as_fil = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_InbreedingCoeff=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_InbreedingCoeff=//;
                                                @AS_Inb = split (/,/, $as);
                                                $as_inb = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_MQ=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_MQ=//;
                                                @AS_MQ = split (/,/, $as);
                                                $as_mq = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_MQRankSum=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_MQRankSum=//;
                                                @AS_MQR = split (/,/, $as);
                                                $as_mqr = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_QD=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_QD=//;
                                                @AS_QD = split (/,/, $as);
                                                $as_qd = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_ReadPosRankSum=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_ReadPosRankSum=//;
                                                @AS_Rea = split (/,/, $as);
                                                $as_rea = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_SOR=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_SOR=//;
                                                @AS_SOR = split (/,/, $as);
                                                $as_sor = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_VQSLOD=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_VQSLOD=//;
                                                @AS_VQS = split (/,/, $as);
                                                $as_vqs = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_culprit=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_culprit=//;
                                                @AS_cul = split (/,/, $as);
                                                $as_cul = $count_info;
                                        }
                    			elsif ($info_data =~ m/^MLEAC=/) {
						my $mleac = $info_data;
						$mleac =~ s/^MLEAC=//;
						@MLEAC = split (/,/, $mleac);
						$mleac_col = $count_info;
		    			}
					elsif ($info_data =~ m/^MLEAF=/) {
						my $mleaf = $info_data;
						$mleaf =~ s/^MLEAF=//;
						@MLEAF = split (/,/, $mleaf);
						$mleaf_col = $count_info;
		    			}

		    			$count_info += 1;
				} # foreach my $info_data (@info)

				my $count_alt = 1;
				foreach my $altern (@alt) {
					my @dados1 = @dados;
		    			my @info1 = @info;
		    			$info1[$ac_col] = "AC=$AC[$count_alt - 1]";
		    			$info1[$af_col] = "AF=$AF[$count_alt - 1]";
					$info1[$as_bas] = "AS_BaseQRankSum=$AS_Bas[$count_alt - 1]";
					$info1[$as_fs] = "AS_FS=$AS_FS[$count_alt - 1]";
					$info1[$as_fil] = "AS_FilterStatus=$AS_Fil[$count_alt - 1]";
					$info1[$as_inb] = "AS_InbreedingCoeff=$AS_Inb[$count_alt - 1]";
					$info1[$as_mq] = "AS_MQ=$AS_MQ[$count_alt - 1]";
					$info1[$as_mqr] = "AS_MQRankSum=$AS_MQR[$count_alt - 1]";
					$info1[$as_qd] = "AS_QD=$AS_QD[$count_alt - 1]";
					$info1[$as_rea] = "AS_ReadPosRankSum=$AS_Rea[$count_alt - 1]";
					$info1[$as_sor] = "AS_SOR=$AS_SOR[$count_alt - 1]";
					$info1[$as_vqs] = "AS_VQSLOD=$AS_VQS[$count_alt - 1]";
					$info1[$as_cul] = "AS_culprit=$AS_cul[$count_alt - 1]";
		    			$info1[$mleac_col] = "MLEAC=$MLEAC[$count_alt - 1]";
		    			$info1[$mleaf_col] = "MLEAF=$MLEAF[$count_alt - 1]";

		    			my $new_info = join (";", @info1);
		    			$dados1[$inf] = $new_info;
		    			$dados1[$alt] = $alt[$count_alt - 1];
		    			my $s_num = $for + 1;

		    			my ($PL1, $PL2);
		    			my ($SAC1, $SAC2);
		    			if ($count_alt == 1) {
						$PL1 = 1;
						$PL2 = 2;
						$SAC1 = 2;
						$SAC2 = 3;
		    			}

		    			elsif ($count_alt == 2) {
                				$PL1 = 3;
                				$PL2 = 5;
						$SAC1 = 4;
						$SAC2 = 5;
		    			}
 
			    		elsif ($count_alt == 3) {
						$PL1 = 6;
						$PL2 = 9;
						$SAC1 = 6;
						$SAC2 = 7;
		    			}

			    		elsif ($count_alt == 4) {
						$PL1 = 10;
						$PL2 = 14;
						$SAC1 = 8;
						$SAC2 = 9;
		    			}

			    		elsif ($count_alt == 5) {
						$PL1 = 15;
						$PL2 = 20;
						$SAC1 = 10;
						$SAC2 = 11;
			    		}	

			    		while ($s_num <= $#dados) {

						if ($dados1[$s_num] =~ m/\.\/\./) {
			    				$dados1[$s_num] = "./.";
						} # if ($dados1[$s_num] =~ m/\.\/\./)
						else {
				    			my @sample = split (/:/, $dados1[$s_num]);
				    			my @GT = split (/\//, $sample[0]);
				    			my @sample_AD = split (/,/, $sample[$AD_col]);
			    				my $AD = "$sample_AD[0],$sample_AD[$count_alt]";
			    				my @sample_PL = split (/,/, $sample[$PL_col]);
			    				my $PL = "$sample_PL[0],$sample_PL[$PL1],$sample_PL[$PL2]";
				    			my $SAC;

				    			if ($sample[$SAC_col] =~ m/[0-9,]+/g) {
								my @sample_SAC = split (/,/, $sample[$SAC_col]);
								$SAC = "$sample_SAC[0],$sample_SAC[1],$sample_SAC[$SAC1],$sample_SAC[$SAC2]";
			    				}
				    			else {
								$SAC = ".";
				    			}

				    			my $DP_s = 0;
			    				foreach my $ADsum (@sample_AD) {
								$DP_s += $ADsum;
			    				}

				    			if ($GT[0] == 0) {

								if ($GT[1] == $count_alt) {
					    				$dados1[$s_num] = "0/1:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
								elsif ($GT[1] == 0) {
				    					$dados1[$s_num] = "0/0:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
								else {
					    				$dados1[$s_num] = "0/.:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}			
			    				} #if ($GT[0] == 0)
				    			elsif (($GT[0] != 0) && ($GT[0] == $GT[1])) {

								if ($GT[0] == $count_alt) {
					    				$dados1[$s_num] = "1/1:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
								else {
				    					$dados1[$s_num] = "./.:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
			    				} # elsif (($GT[0] != 0) && ($GT[0] == $GT[1]))
				    			elsif (($GT[0] != 0) && ($GT[0] != $GT[1])) {

								if ($GT[0] == $count_alt) {
					    				$dados1[$s_num] = "./1:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
								elsif ($GT[1] == $count_alt) {
				    					$dados1[$s_num] = "./1:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
								else {
				    					$dados1[$s_num] = "./.:$AD:$DP_s:$sample[$GQ_col]:$PL:$SAC";
								}
			    				} # elsif (($GT[0] != 0) && ($GT[0] != $GT[1]))
						} # else

						$s_num += 1;
		    			} # while ($s_num <= $total_col)

			    		my $new_line = join ("\t", @dados1);
			    		print OUT "$new_line\n";
			    		$totalLinePos += 1;
		    			$count_alt += 1;
				} # foreach my $altern (@alt)
    

		    	} #if ($num_alt > 1)
		    	else {
				my $s_num = $for + 1;

				while ($s_num <= $#dados) {
	
				    	if ($dados[$s_num] =~ m/\.\/\./) {
						$dados[$s_num] = "./.";
						$s_num += 1;
		    			} # if ($dados[$s_num] =~ m/\.\/\./) 
		    			else {
						my @sample = split (/:/, $dados[$s_num]);
						my @sample_AD = split (/,/, $sample[$AD_col]);

						my $SAC;
						if ($sample[$SAC_col] =~ m/[0-9,]+/g) {
				    			$SAC = $sample[$SAC_col];
						}
						else {
			    				$SAC = ".";
						}
			
						my $DP_s = 0;
						foreach my $AD_sum (@sample_AD) {
				    			$DP_s += $AD_sum;
						}

						$dados[$s_num] = "$sample[0]:$sample[$AD_col]:$DP_s:$sample[$GQ_col]:$sample[$PL_col]:$SAC";
						$s_num += 1;
		    			} # else
				} # while ($s_num <= $total_col)

				my $new_line = join ("\t", @dados);
				print OUT "$new_line\n";
				$totalLinePos += 1;
				$totalLineSingle += 1;
	    		} # else
		} # if ($exist_PGT == 0)

		elsif ($exist_PGT == 1) {
			$dados[$for] = "GT:AD:DP:GQ:PL:PGT:PID:SAC";
			if ($num_alt > 1) {
				$totalLineSplit += 1;
				my @info = split(/;/, $dados[$inf]);
	
				my (@AC, @AF);
                                my (@AS_Bas, @AS_FS, @AS_Fil, @AS_Inb, @AS_MQ, @AS_MQR, @AS_QD, @AS_Rea, @AS_SOR, @AS_VQS, @AS_cul);
                                my (@MLEAC, @MLEAF);
                                my ($ac_col, $af_col);
                                my ($as_bas, $as_fs, $as_fil, $as_inb, $as_mq, $as_mqr, $as_qd, $as_rea, $as_sor, $as_vqs, $as_cul);
                                my ($mleac_col, $mleaf_col);
	
				my $count_info = 0;
				foreach my $info_data (@info) {
					if ($info_data =~ m/^AC=/) {
						my $ac = $info_data;
						$ac =~ s/^AC=//;
						@AC = split (/,/, $ac);
						$ac_col = $count_info;
		    			}
				    	elsif ($info_data =~ m/^AF/) {
						my $af = $info_data;
						$af =~ s/^AF=//;
						@AF = split (/,/, $af);
						$af_col = $count_info;
		    			}
					elsif ($info_data =~ m/^AS_BaseQRankSum=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_BaseQRankSum=//;
                                                @AS_Bas = split (/,/, $as);
                                                $as_bas = $count_info;
                                        }               
                                        elsif ($info_data =~ m/^AS_FS=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_FS=//;
                                                @AS_FS = split (/,/, $as);
                                                $as_fs = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_FilterStatus=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_FilterStatus=//;
                                                @AS_Fil = split (/,/, $as);
                                                $as_fil = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_InbreedingCoeff=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_InbreedingCoeff=//;
                                                @AS_Inb = split (/,/, $as);
                                                $as_inb = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_MQ=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_MQ=//;
                                                @AS_MQ = split (/,/, $as);
                                                $as_mq = $count_info;
                                        }
					elsif ($info_data =~ m/^AS_MQRankSum=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_MQRankSum=//;
                                                @AS_MQR = split (/,/, $as);
                                                $as_mqr = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_QD=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_QD=//;
                                                @AS_QD = split (/,/, $as);
                                                $as_qd = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_ReadPosRankSum=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_ReadPosRankSum=//;
                                                @AS_Rea = split (/,/, $as);
                                                $as_rea = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_SOR=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_SOR=//;
                                                @AS_SOR = split (/,/, $as);
                                                $as_sor = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_VQSLOD=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_VQSLOD=//;
                                                @AS_VQS = split (/,/, $as);
                                                $as_vqs = $count_info;
                                        }
                                        elsif ($info_data =~ m/^AS_culprit=/) {
                                                my $as = $info_data;
                                                $as =~ s/^AS_culprit=//;
                                                @AS_cul = split (/,/, $as);
                                                $as_cul = $count_info;
                                        }
					elsif ($info_data =~ m/^MLEAC=/) {
						my $mleac = $info_data;
						$mleac =~ s/^MLEAC=//;
						@MLEAC = split (/,/, $mleac);
						$mleac_col = $count_info;
					}
					elsif ($info_data =~ m/^MLEAF=/) {
						my $mleaf = $info_data;
						$mleaf =~ s/^MLEAF=//;
						@MLEAF = split (/,/, $mleaf);
						$mleaf_col = $count_info;
		    			}

					$count_info += 1;
				} # foreach my $info_data (@info)
		
				my $count_alt = 1;
				foreach my $altern (@alt) {
					my @dados1 = @dados;
					my @info1 = @info;
					$info1[$ac_col] = "AC=$AC[$count_alt - 1]";
					$info1[$af_col] = "AF=$AF[$count_alt - 1]";
					$info1[$as_bas] = "AS_BaseQRankSum=$AS_Bas[$count_alt - 1]";
                                        $info1[$as_fs] = "AS_FS=$AS_FS[$count_alt - 1]";
                                        $info1[$as_fil] = "AS_FilterStatus=$AS_Fil[$count_alt - 1]";
                                        $info1[$as_inb] = "AS_InbreedingCoeff=$AS_Inb[$count_alt - 1]";
                                        $info1[$as_mq] = "AS_MQ=$AS_MQ[$count_alt - 1]";
                                        $info1[$as_mqr] = "AS_MQRankSum=$AS_MQR[$count_alt - 1]";
                                        $info1[$as_qd] = "AS_QD=$AS_QD[$count_alt - 1]";
                                        $info1[$as_rea] = "AS_ReadPosRankSum=$AS_Rea[$count_alt - 1]";
                                        $info1[$as_sor] = "AS_SOR=$AS_SOR[$count_alt - 1]";
                                        $info1[$as_vqs] = "AS_VQSLOD=$AS_VQS[$count_alt - 1]";
                                        $info1[$as_cul] = "AS_culprit=$AS_cul[$count_alt - 1]";
					$info1[$mleac_col] = "MLEAC=$MLEAC[$count_alt - 1]";
					$info1[$mleaf_col] = "MLEAF=$MLEAF[$count_alt - 1]";

					my $new_info= join (";", @info1);
					$dados1[$inf] = $new_info;
					$dados1[$alt] = $alt[$count_alt - 1];
					my $s_num = $for + 1;

					my ($PL1, $PL2);
					my ($SAC1, $SAC2);
					if ($count_alt == 1) {
						$PL1 = 1;
						$PL2 = 2;
						$SAC1 = 2;
						$SAC2 = 3;
					} # if ($count_alt == 1)
					elsif ($count_alt == 2) {
						$PL1 = 3;
						$PL2 = 5;
						$SAC1 = 4;
						$SAC2 = 5;
		    			}
					elsif ($count_alt == 3) {
						$PL1 = 6;
						$PL2 = 9;
						$SAC1 = 6;
						$SAC2 = 7;
		    			}
					elsif ($count_alt == 4) {
						$PL1 = 10;
						$PL2 = 14;
						$SAC1 = 8;
						$SAC2 = 9;
		    			}
					elsif ($count_alt == 5) {
						$PL1 = 15;
						$PL2 = 20;
						$SAC1 = 10;
						$SAC2 = 11;
		    			}

					while ($s_num <= $#dados) {
	
						if ($dados1[$s_num] =~ m/\.\/\./) {
							$dados1[$s_num] = "./.";
						} # if ($dados1[$s_num] =~ m/\.\/\./)
						else {
							my @sample = split (/:/, $dados1[$s_num]);
							my @GT = split (/\//, $sample[0]);
							my @sample_AD = split (/,/, $sample[$AD_col]);
							my $AD = "$sample_AD[0],$sample_AD[$count_alt]";
							my @sample_PL = split (/,/, $sample[$PL_col]);
							my $PL = "$sample_PL[0],$sample_PL[$PL1],$sample_PL[$PL2]";
							my $SAC;
							if ($sample[$SAC_col] =~ m/[0-9,]+/g) {
								my @sample_SAC = split (/,/, $sample[$SAC_col]);
								$SAC = "$sample_SAC[0],$sample_SAC[1],$sample_SAC[$SAC1],$sample_SAC[$SAC2]";
			    				}
			    				else {
								$SAC= ".";
			    				}

							my $DP_s = 0;
							foreach my $ADsum (@sample_AD) {
								$DP_s += $ADsum;
			    				}

							if ($GT[0] == 0) {
								if ($GT[1] == $count_alt) {
									$dados1[$s_num] = "0/1:$AD:$DP_s:$sample[$GQ_col]:$PL:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
								} # if ($GT[1] == $count_alt) 
								elsif ($GT[1] == 0) {
									$dados1[$s_num] = "0/0:$AD:$DP_s:$sample[$GQ_col]:$PL:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
								}
								else {
									$dados1[$s_num] = "0/.:$AD:$DP_s:$sample[$GQ_col]:$PL:.:.:$SAC";
								}
			    				} # if ($GT[0] == 0)
							elsif (($GT[0] != 0) && ($GT[0] == $GT[1])) {
								if ($GT[0] == $count_alt) {
									$dados1[$s_num] = "1/1:$AD:$DP_s:$sample[$GQ_col]:$PL:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
								}
								else {
									$dados1[$s_num] = "./.:$AD:$DP_s:$sample[$GQ_col]:$PL:.:.:$SAC";
								}
			    				} # elsif (($GT[0] != 0) && ($GT[0] == $GT[1]))
			    				elsif (($GT[0] != 0) && ($GT[0] != $GT[1])) {
								if ($GT[0] == $count_alt) {
									$dados1[$s_num] = "./1:$AD:$DP_s:$sample[$GQ_col]:$PL:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
								}
								elsif ($GT[1] == $count_alt) {
									$dados1[$s_num] = "./1:$AD:$DP_s:$sample[$GQ_col]:$PL:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
								}
								else {
									$dados1[$s_num] = "./.:$AD:$DP_s:$sample[$GQ_col]:$PL:.:.:$SAC";
								} 
			    				} # elsif (($GT[0] != 0) && ($GT[0] != $GT[1]))
						} # else

						$s_num += 1;
		    			} # while ($s_num <= $total_col)
		
					my $new_line = join ("\t", @dados1);
					print OUT "$new_line\n";
					$totalLinePos += 1;
					$count_alt += 1;

				} #foreach my $altern (@alt)

	    		} #if ($num_alt > 1)	

			else {
				my $s_num = $for + 1;
				while ($s_num <= $#dados) {

					if ($dados[$s_num] =~ m/\.\/\./) {
						$dados[$s_num] = "./.";
		    			} # if ($dados[$s_num] =~ m/\.\/\./)
					else {
						my @sample = split (/:/, $dados[$s_num]);
						my @sample_AD = split (/,/, $sample[$AD_col]);

						my $SAC;
						if ($sample[$SAC_col] =~ m/[0-9,]+/g) {
							$SAC = $sample[$SAC_col];
						} # if ($sample[$SAC_col] =~ m/[0-9,]+/g)
						else {
							$SAC = ".";
						}

						my $DP_s = 0;
						foreach my $AD_sum (@sample_AD) {
							$DP_s += $AD_sum;
						}

						$dados[$s_num] = "$sample[0]:$sample[$AD_col]:$DP_s:$sample[$GQ_col]:$sample[$PL_col]:$sample[$PGT_col]:$sample[$PID_col]:$SAC";
		    			} # else

					$s_num += 1;
				} # while ($s_num <= $total_col) 
		
				my $new_line = join ("\t", @dados);
				print OUT "$new_line\n";
				$totalLinePos += 1;
				$totalLineSingle += 1;
	    		} # else 
		} #elsif ($exist_PGT == 1)
	} #elsif ($line =~ m/^chr/)
}


close (IN);
close (OUT);

print "Lines\tTotal input/Total single/Total split/Total output: $totalLinePre/$totalLineSingle/$totalLineSplit/$totalLinePos\n";

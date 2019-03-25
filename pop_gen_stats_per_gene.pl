## Only includes non-polarized stats


use strict;
use warnings;
use POSIX;


my $min_call_rate = 0.5;
my @populations = (qw/Pop1 Pop2/);

open OUT, ">Pop_gen_stats_per_gene.txt" or die;
open GENE, "Pfal_CDS_ns_syn_ffd_count_r.31.txt" or die;
while (<GENE>) {
  chomp;
  my @gene_line = split("\t", $_);

### Make output file header
  if ($. == 1) {
    print OUT join("\t", @gene_line[0..14]);
    foreach my $pop (@populations) {
      print OUT "\t$pop"."_S\t$pop"."_pi\t$pop"."_ThetaW\t$pop"."_TD";
    }
    foreach my $index1 (0..($#populations-1)) {
      my $pop1 = $populations[$index1];
      foreach my $index2 (($index1+1)..$#populations) {
	my $pop2 = $populations[$index2];
	print OUT "\t$pop1"."-$pop2"."_shared_alleles\tFST_$pop1"."-$pop2"."\tTr_FST_$pop1"."-$pop2";
      }
    }
    print OUT "\n";
    next;
  }
  my $gene = $gene_line[0];
  foreach my $numeral (1..9) {
    $gene =~ s/\.$numeral//g;
  }
  my $chrom = $gene_line[1];
  my $site_number = $gene_line[3]; ## Number of NS sites per trans
  my $length = $site_number;
  my $wl = $site_number;
  if ($wl == 0) {
    print "Gene $gene (line $.) has length of 0\n";
    next;
  }
  my $gene_info = join("\t", @gene_line[0..14]);

  my @fst_parts_columns;
  my @pi_columns;
  my @ThetaW_columns;
  my @Sample_Coverage_columns;
  my @Alt_Count_columns;

  open SNPS, "$gene"."_pop_gen_stats_per_nt.txt" or (print "No gene variant file for gene $gene\n" and next);
  
  my %fst_a_sum = my %fst_b_sum = my %fst_c_sum = my %snp_count = my %pi_sum = my %ThetaW_sum = my %a1_sum = my %a2_sum = my %b1_sum  = my %b2_sum  = my %c2_sum  = my %ThetaH_sum = my %ThetaL_sum = my %FayWuN_sum = my %bn1_sum = my %vd_sum = my %ud_sum = my %ne_sum = my %KST_sum = my %Snn_sum = my %Sample_Coverage = ();
  my $total_snps = 0;
  my $snps_remain = 0;
  while (<SNPS>) {
    chomp;
    my @dl2 = split("\t", $_);
    if ($.==1) {
      foreach my $col (0..$#dl2) {
	if ($dl2[$col] =~ /_PARTS/) {
	  @fst_parts_columns = (@fst_parts_columns, $col);
	}
	if ($dl2[$col] =~ /_pi/) {
	  @pi_columns = (@pi_columns, $col);
	}
	if ($dl2[$col] =~ /_WatTheta/) {
	  @ThetaW_columns = (@ThetaW_columns, $col);
	}
	if ($dl2[$col] =~ /CALLS_/) {
	  @Sample_Coverage_columns = (@Sample_Coverage_columns, $col);
	}
	if ($dl2[$col] =~ /ALT_COUNT_/) {
	  @Alt_Count_columns = (@Alt_Count_columns, $col);
	}
      }
#      print OUT "$gene_info";
      next;
    }

    if ($dl2[4] ne "PASS") {
      next;
    }

    $snps_remain = 1;
    
## Count number of segregating sites
    ++$total_snps;
    
## Make hash of FST parts
    foreach my $pair (@fst_parts_columns) {
      if ($dl2[$pair] eq "NA") {
	next;
      }
      if (!exists $snp_count{$pair}) {
	$fst_a_sum{$pair} = 0;
	$fst_b_sum{$pair} = 0;
	$fst_c_sum{$pair} = 0;
	$snp_count{$pair} = 0;
      }
      my @fst_parts = split(";",$dl2[$pair]);
      $fst_a_sum{$pair} = ($fst_a_sum{$pair} + $fst_parts[0]);
      $fst_b_sum{$pair} = ($fst_b_sum{$pair} + $fst_parts[1]);
      $fst_c_sum{$pair} = ($fst_c_sum{$pair} + $fst_parts[2]);
      ++$snp_count{$pair};
    }
    
## Make Hash of pi parts
    foreach my $pi (@pi_columns) {
      if (!exists $pi_sum{$pi}){
	$pi_sum{$pi} = 0;
	$snp_count{$pi} = 0;
      }
      if ($dl2[$pi] ne "NA") {
	$pi_sum{$pi} = ($pi_sum{$pi} + $dl2[$pi]);
	if ($dl2[$pi] > 0) {
	  ++$snp_count{$pi};
	}
      }
    }
    
## Make Hash of ThetaW parts
    foreach my $index (0..$#ThetaW_columns) {
    my $ThetaW = $ThetaW_columns[$index];
      if (!exists $ThetaW_sum{$ThetaW}) {
	$ThetaW_sum{$ThetaW} = 0;
	$snp_count{$ThetaW} = 0;
      }
      my @theta_parts = split(";",$dl2[$ThetaW]);
      if ($theta_parts[0] ne "NA") {
	$ThetaW_sum{$ThetaW} = ($ThetaW_sum{$ThetaW} + $theta_parts[0]);
	if ($theta_parts[0] > 0) {
	  my $count_col = $Alt_Count_columns[$index];
	  my @alt_counts = split(",", $dl2[$count_col]);
	  my $alt_allele_count = 0;
	  foreach my $allele (@alt_counts) {
	    if ($allele > 0) {
	      ++$alt_allele_count;
	    }
	  }
	  $snp_count{$ThetaW} = $snp_count{$ThetaW} + $alt_allele_count;
	}
      }
    }

## Count cumulative sample coverage across the sites
    foreach my $index (0..$#pi_columns) {
      my $pop = $Sample_Coverage_columns[$index];
      if (!exists $Sample_Coverage{$pop}) {
	$Sample_Coverage{$pop} = 0;
      }
      if ($dl2[$pi_columns[$index]] ne "NA") {
	if ($dl2[$pi_columns[$index]] > 0) {
	  $Sample_Coverage{$pop} = $Sample_Coverage{$pop} + $dl2[$pop];
	}
      }
    }
  }
  close SNPS;

  if ($snps_remain == 0) { ## all snps in file were filtered out
    next;
  } else {
    print OUT "$gene_info";
  }
   
  my $avg_pi = 0;
  my $ThetaW = 0;
  my $TD = "NA";
  foreach my $index (0..$#pi_columns) { ## cycle through individual population statistics
## Calculate Average pi
    my $pi_col = $pi_columns[$index];
    my $Watterson_col = $ThetaW_columns[$index];
    my $coverage_col = $Sample_Coverage_columns[$index];
    if (!exists $pi_sum{$pi_col}) {
      $pi_sum{$pi_col} = 0;
    }
    my $avg_pi = ($pi_sum{$pi_col}/$wl);
    
    ## Calculate Watterson's Theta and Tajima's D
    if (!exists $snp_count{$pi_col}) {
      $snp_count{$pi_col} = 0;
    }
    if (!exists $snp_count{$Watterson_col}) {
      $snp_count{$Watterson_col} = 0;
    }
    if ($snp_count{$pi_col} > 0) {
      $ThetaW = $ThetaW_sum{$Watterson_col}/$wl;
      my $sample = ceil($Sample_Coverage{$coverage_col}/$snp_count{$Watterson_col}); ### calculates average sample size across all the sites
      my $a1 =my $a2 = 0;
      foreach my $n_TD (1..($sample-1)) {
	$a1 = ($a1 + (1/$n_TD));
	$a2 = ($a2 + (1/($n_TD*$n_TD)));
      }
      my $b1 =($sample+1)/(3*($sample-1));
      my $b2 = ((2*(($sample*$sample)+$sample+3))/(9*$sample*($sample-1)));
      my $c1 = $b1-(1/$a1);
      my $c2 = ($b2-(($sample+2)/($a1*$sample))+($a2/($a1*$a1)));
      my $e1 = $c1/$a1;
      my $e2 = $c2/(($a1*$a1)+$a2);
      my $S = ($snp_count{$pi_col}/$wl);
      my $TD_var = (sqrt(($e1*$snp_count{$pi_col})+($e2*$snp_count{$pi_col}*($snp_count{$pi_col}-1))));
      if ($TD_var != 0) {
	$TD = (($pi_sum{$pi_col}-($snp_count{$Watterson_col}/$a1))/$TD_var);
      }
   }

## Print individual population stats
    print OUT "\t$snp_count{$pi_col}\t$avg_pi\t$ThetaW\t$TD";
  }

## Average FST
  foreach my $index (0.. $#fst_parts_columns) {
    my $avg_fst ="NA";
    my $Tr_FST = "NA";
    my $pair_col = $fst_parts_columns[$index];
    if (!exists $snp_count{$pair_col}) {
      $snp_count{$pair_col} = 0;
    }
    if ($snp_count{$pair_col} > 0) {
      if (($fst_a_sum{$pair_col} + $fst_b_sum{$pair_col} + $fst_c_sum{$pair_col}) > 0) {
	$avg_fst = ($fst_a_sum{$pair_col}/($fst_a_sum{$pair_col} + $fst_b_sum{$pair_col} + $fst_c_sum{$pair_col}));
	$Tr_FST = $avg_fst;
	if ($avg_fst < 0) {
	  $Tr_FST = 0;
	}
      }
    }
    print OUT "\t$snp_count{$pair_col}\t$avg_fst\t$Tr_FST";
  }
  print OUT "\n";
}
close GENE;
close OUT;	

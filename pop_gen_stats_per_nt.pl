## Multiallelic
## No polarization
## No minimum sample size limit

use warnings;
use strict;

use List::Util qw(sum);


my $min_call_rate = 0.5;

my @populations = (qw/Pop1 Pop2/);


###### List of gene names (or file prefix) used for each individual input file. One ID per line.
my @gene;
open LIST, "gene_IDs.txt" or die;
while (<LIST>) {
  chomp;
  @gene = (@gene, $_);
}
close LIST;


foreach my $gene (@gene) {
my $input = "$gene"."_summary_filtered_samples.txt";
my $output = "$gene"."_pop_gen_stats_per_nt.txt";

  
  open SNPS, $input or (print "No summary snp file for gene $gene\n" and next);
  open OUT, ">$output";
  my %POS;
  my %TOTAL;
  my %ALT_COL;
  my %REF_COL;
  my %CALLED_COL;
  my %PERC_CALLED_COL;
  my %MISSING_COL;
  while (<SNPS>) {
    chomp;
    my (@dl) = split("\t",$_);
    if ($. == 1){
      print OUT "$_";
      foreach my $pop (@populations) {
	print OUT "\t$pop"."_pi\t$pop"."_WatTheta";
      }
      print OUT "\tTotal_allele_counts";
      foreach my $index1 (0..($#populations-1)) {
	my $pop1 = $populations[$index1];
	foreach my $index2 (($index1+1)..$#populations) {
	  my $pop2 = $populations[$index2];
	  print OUT "\t$pop1"."-$pop2"."_shared_alleles\tFST_$pop1"."-$pop2"."\tTr_FST_$pop1"."-$pop2\t$pop1"."-$pop2"."_PARTS";
	}
      }
      print OUT "\n";
      foreach my $col (0..$#dl) {
	if ($dl[$col] =~ /ALT_COUNT/) {
	  my ($other, $pop) = split("_COUNT_",$dl[$col]);
	  $ALT_COL{$pop} = $col;
	}
	if ($dl[$col] =~ /REF_COUNT/) {
	  my ($other, $pop) = split("_COUNT_",$dl[$col]);
	  $REF_COL{$pop} = $col;
	}
	if ($dl[$col] =~ /CALLS/) {
	  my ($other, $pop) = split("_",$dl[$col]);
	  $CALLED_COL{$pop} = $col;
	}
	if ($dl[$col] =~ /MISSING/) {
	  my ($other, $pop) = split("_",$dl[$col]);
	  $MISSING_COL{$pop} = $col;
	}
	if ($dl[$col] =~ /PERCENT_CALLED/) {
	  my ($other, $pop) = split("_CALLED_",$dl[$col]);
	  $PERC_CALLED_COL{$pop} = $col;
	}
      }
      next;
    }
    
    my %COUNTS;
    my $site = $dl[1];
    my $Ref_allele = $dl[1];
    my $Alt_allele = $dl[3];
    my @alleles = split(",",$Alt_allele);
    my $allele_count = ($#alleles+2);
    my @total_count = ();
    foreach (0..($allele_count-1)) {
      $total_count[$_] = 0;
    }
    print OUT "$_";
  
## Calculate unpolarized diversity statistics
    my %CALLED_SITES;
    foreach my $pop (@populations) {
      if ($dl[$PERC_CALLED_COL{$pop}] < $min_call_rate) {
	print OUT "\tNA\tNA";
	$COUNTS{$pop} = "NA";
      } else {
	my @alt_alleles;
	my $alt_count = 0;
	if ($allele_count > 2) {
	  @alt_alleles = split(",",$dl[$ALT_COL{$pop}]);
	  foreach my $allele (@alt_alleles) {
	    $alt_count += $allele;
	  }
	} else {
	  @alt_alleles = ($dl[$ALT_COL{$pop}]);
	  $alt_count = $dl[$ALT_COL{$pop}];
	}
	my $ref_allele = $dl[$REF_COL{$pop}];
	my $called_sites = $ref_allele+$alt_count;
	$CALLED_SITES{$pop} = $called_sites;
	my @all_alleles = ($ref_allele,@alt_alleles);
	foreach my $index (0..$#all_alleles) {
	  $total_count[$index] = $total_count[$index] + $all_alleles[$index];
	}
      
## pi
	my $pi;
	$COUNTS{$pop} = "$ref_allele,$dl[$ALT_COL{$pop}]"; # comma separated list of allele frequencies within the population
	if (($alt_count == 0) || (($ref_allele==0) && ($allele_count==2))) {
	  $pi = "NA";
	} else {
	  my $hetero_sum = 1;
	  foreach my $allele (@all_alleles) {
	    my $freq = $allele/$called_sites;
	    $hetero_sum = $hetero_sum - ($freq*$freq);
	  }
	  $pi = ($called_sites/($called_sites-1))*($hetero_sum);
	}
      
## Watterson's Theta and Tajima's D
	my $a1_sum = my $a2_sum = my $b1_sum = my $b2_sum = my $c2_sum = 0;
	my $WatTheta_sum  = 0;
	if ($pi ne "NA") {
	  foreach my $allele (@alt_alleles) {
	    my $a1_allele = my $a2_allele = my $b1_allele = my $b2_allele = my $c2_allele = my $WatTheta_allele = 0;
	    if ($allele > 0) {
	      foreach my $index (1..($called_sites-1)) {
		$a1_allele = $a1_allele+(1/$index);
		$a2_allele = $a2_allele+(1/($index*$index));
	      }
	      $WatTheta_allele = (1/$a1_allele);
	    }
	    $WatTheta_sum = $WatTheta_sum + $WatTheta_allele;
	  }
	}
	print OUT "\t$pi\t$WatTheta_sum";

      }
    }
      
### Calculate pairwise FST
    my $total_counts = join(";",@total_count);
    print OUT "\t$total_counts";
    foreach my $first (0..($#populations-1)){
      my $first_pop = $populations[$first];
      foreach my $second (($first + 1)..$#populations) {
	my $second_pop = $populations[$second];
	&FST($allele_count,$COUNTS{$first_pop},$COUNTS{$second_pop},$CALLED_SITES{$first_pop},$CALLED_SITES{$second_pop},$total_counts);

      }
    }
    print OUT "\n";
  }
close OUT;
close SNPS;
}

###### Multiallelic FST, Uses Robertson and Hill (1984) approach to weighting, which Weir and Cockerham (1984) show to be approporiate for populations with low FST
sub FST {
  if (($_[1] =~ /NA/) || ($_[2] =~ /NA/)) {
    print OUT "\tNA\tNA\tNA\tNA";
    return;
  }
  my $r = 2;
  
  my $allele_count = $_[0]; ########### Check whether all alleles are segregating in the populations
  my $shared_allele_count = $allele_count;
  foreach my $rep (0..($allele_count-1)) { # Calculate ThetaWC for each allele
    my @count1 = split(",",$_[1]);  ##### Check whether all alleles are segregating in the populations
    my @count2 = split(",",$_[2]);
    if (($count1[$rep] == 0) && ($count2[$rep]==0)) {
      $shared_allele_count = $shared_allele_count - 1;
    }
  }

  if ($shared_allele_count == 1) {
    print OUT "\t1\tNA\tNA\tNA";
    return;
  }
  
  my @anc_count = split(";",$_[5]);
  my $avg_n = (($_[3] + $_[4])/2); # average of the two population's sample sizes
  my $weighted_fst = 0;
  my $fst_a_weighted = 0;
  my $fst_b_weighted = 0;
  my $fst_c_weighted = 0;
  foreach my $rep (0..($allele_count-1)) { # Calculate ThetaWC for each allele
    my @count1 = split(",",$_[1]);  ##### only include alleles that are segregating in all populations
    my @count2 = split(",",$_[2]);
    if (($count1[$rep] == 0) && ($count2[$rep]==0)) {
      next;
    }
    my @n_values = my @nc_values = my @p_freq_values = my @p_pop_values = my @h_values = ();
    foreach my $pop (1,2) {
      my @count = split(",",$_[$pop]);
      my $p_count = $count[$rep]; # number of target alleles in the population
      my $pop_n = $_[($pop +2)]; # total sample size of the population
      my $p_pop_freq = ($p_count/$pop_n); # frequency of the allele in the populaiton
      my $p_pop = (($pop_n*$p_pop_freq)/($r*$avg_n));
      my $nc_pop = (($pop_n*$pop_n)/($r*$avg_n));
      my $h_pop = (($pop_n*(2*($p_pop_freq)*(1-$p_pop_freq)))/(($r*$avg_n)));
      push(@nc_values, $nc_pop);
      push(@p_freq_values, $p_pop_freq);
      push(@p_pop_values, $p_pop);
      push(@h_values, $h_pop);
      push(@n_values, $pop_n);
    }    
    my $sum_nc_pops = ($nc_values[0] + $nc_values[1]);
    my $p = ($p_pop_values[0] + $p_pop_values[1]);
    my $h = ($h_values[0] + $h_values[1]);
    my $total_n = ($n_values[0] + $n_values[1]);
    my $avg_p_freq = (($p_freq_values[0]+$p_freq_values[1])/$r);
    my $avg_n = ($total_n/$r);
    
    my $s2 = ((($n_values[0]*($p_freq_values[0]-$avg_p_freq)*($p_freq_values[0]-$avg_p_freq))/(($r-1)*$avg_n)) + (($n_values[1]*($p_freq_values[1]-$avg_p_freq)*($p_freq_values[1]-$p))/(($r-1)*$avg_n)));
    my $nc = ((($r*$avg_n)-$sum_nc_pops)/($r-1));
##########
    if ($nc == 0) {
      print OUT "\t$allele_count\tNA\tNA\tNA";
      return;
    } else {
      ############
      my $combined_pop_count = $_[3]+$_[4];
      my $weight = (1 - ($anc_count[$rep]/$combined_pop_count))/($allele_count - 1);
      my $fst_a = (($avg_n/$nc)*($s2-(1/($avg_n-1))*($p*(1-$p)-((($r-1)/$r)*$s2)-($h/4))));
      my $fst_b = (($avg_n/($avg_n-1))*($p*(1-$p)-((($r-1)/$r)*$s2)-((((2*$avg_n)-1)/(4*$avg_n))*$h)));
      my $fst_c = ($h/2);
      my $raw_fst = $fst_a/($fst_a+$fst_b+$fst_c);
      $weighted_fst = $weighted_fst + ($weight * ($fst_a/($fst_a+$fst_b+$fst_c)));
      $fst_a_weighted = $fst_a_weighted + ($weight * $fst_a);
      $fst_b_weighted = $fst_b_weighted + ($weight * $fst_b);
      $fst_c_weighted = $fst_c_weighted + ($weight * $fst_c);
    }
  }
  my $FST_Tr = $weighted_fst;
  if ($weighted_fst < 0) {
    $FST_Tr = 0;
  }
  if ($weighted_fst > 1) {
    print "PROBLEM! Found pairwise FST of $weighted_fst for $_[3] and $_[4]\n";
  }
  print OUT "\t$allele_count\t$weighted_fst\t$FST_Tr\t$fst_a_weighted;$fst_b_weighted;$fst_c_weighted";
}

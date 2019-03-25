# Dantzler_et_al_Diversity_Calcs
Scripts and example input files used to calculate diversity statistics for Dantzler et al, 2019

##################
##################Sample Input Files
##################
#SampleGene_summary_filtered_samples.txt
Sample input file for pop_gen_stats_per_nt.pl
Lists variants in genomic region of interest and summarizes allele counts for each population

#Pfal_CDS_ns_syn_ffd_count_r.31.txt
Lists counts of nonsynonymous, synonymous, and four-fold degenerate sites for each annotated P. falciparum transcript (PlasmoDB 3D7 r.31)
Used in pop_gen_stats_per_gene.pl


##################
################## Scripts
##################

#pop_gen_stats_per_nt.pl
Input: SampleGene_summary_filtered_samples.txt, gene_IDs.txt 
Output: SampleGene_pop_gen_stats_per_nt.txt"

#pop_gen_stats_per_gene.pl
Input: SampleGene_pop_gen_stats_per_nt.txt (pop_gen_stats_per_nt.pl output), Pfal_CDS_ns_syn_ffd_count_r.31.txt
Output: Pop_gen_stats_per_gene.txt

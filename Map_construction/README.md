#Files for map construction in Lep-Map 3 and Lep-Anchor. Also includes scripts for converting Lep-Map genotype output for use in R/QTL2.
Full_map_construction_script_May_selfing_final_clean.txt - commands to call genotype probabilities, split linkage groups, order markers, identify highest likelihood marker order, trim bad markers, associate linkage map with draft assembly, export genotypes, process genotypes for use in R/QTL2, and export a fasta file of the pseudomolecule assembly. 
Extract_genotypes.py - approximates raw base calls from the Lep-Map post.call file to assign alleles to grandparents.
Ipomoea_map_LG_position_switcher_Lep_Anchor.py - convert SNP positions from draft assembly position to pseudomolecule position, using .agp file and linkage map. 
full_genome.agp - .agp file relating draft assembly to pseudomolecules based on linkage groups.

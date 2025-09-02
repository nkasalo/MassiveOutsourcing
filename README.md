# Massive Outsourcing of energetically costly amino acids at the origin of animals

The repository contains the sample data and the analysis scripts accompanying the manuscript Kasalo, N., Domazet-Lošo, M., Domazet-Lošo, T. **Massive outsourcing of energetically costly amino acids at the origin of animals**.

All code was written in RStudio version 2025.05.0+496, R version 4.3.2, on macOS Monterey version 12.0.1.
To run the amino acid auxotrophy detection protocol, MMSeqs2 has to be installed (tested on versions 14-7e284 and 15-6f452).


The repository contains two protocols:
1. Combinatorial phenotype selection test
2. Amino acid auxotrophy detection


**Combinatorial phenotype selection test**

The provided example contains the following files:
1. aa_labels.txt - names of amino acids and their three- and one-letter codes
2. disp_cost_table.xlsx - energetic costs of each amino acid and their canonical dispensability
3. listAA20_energy_new_all.txt - energetic costs of each amino acid
4. combinatorial phenotype selection test.R - the R script used to perform the test

The example is set up to calculate the p-value and generate the corresponding graph for the mean opportunity cost of the observed set of 11 EAAs. The code to generate other permutations is provided in the beginning of the file.


**Amino acid auxotrophy detection**

This directory contains the following files:
1. AA_enzymes_eukarya_sample.faa - fasta file of reference enzymes. NOTE: this is a shortened sample; complete data is available in supplementary data of this manuscript
2. aa_labels.txt - names of amino acids and their three- and one-letter codes
3. Biosynthethic_pathways.xlsx - Data on enzymes involved in amino acid biosynthethic pathways
4. DB_2024_final - a directory containing a sample of proteomes in fasta format, named according to the species' taxIDs. The full database is available in Domazet-Lošo, M., Široki, T., Šimičević, K. et al. Macroevolutionary dynamics of gene family gain and loss along multicellular eukaryotic lineages. Nat Commun 15, 2663 (2024). https://doi.org/10.1038/s41467-024-47017-w
5. heatmap.R - the R script used to assign detected enzymes to pathways and draw the heatmap
6. listAA20_energy_new_all.txt - energetic costs of each amino acid
7. mmseq_individual_clustering.R - the R script used to perform sequence search of proteomes against the reference database of enzymes to recover homologous genes
8. names_2023.txt - mapping of species names to taxIDs
9. pathways_all_KEGG_organisms_eukarya.tsv - mapping of genes to amino acid biosynthesis pathway data
10. phylo_2023.xml - the phylogenetic tree in XML format
11. examples

The analysis is performed by first running the mmseq_individual_clustering.R script, which uses MMSeqs2 to identify homologs of genes provided in the reference database in the proteomes of all given species. There are two protocols that can be used: reciprocal best hits or clustering, both of which produce a file that shows the mapping of genes from the given species to genes in the reference database. Both approaches analyze one species at a time, which may take between a few minutes and half an hour, depending on the e-value cutoff, genome size, reference database size, and the chosen protocol.

For reciprocal best hits, the output format is the standard BLAST output, available as mmseq_biosynthesis_reciprocal_final.tsv in the "examples" directory. The clustering output contains only the first two columns which represent the mapping of genes, available as mmseq_biosynthesis_final_08_e-40.tsv in the "examples" directory.

This file is then processed using the heatmap.R script, which for each species assigns each detected enzyme to its corresponding biosynthethic pathway, calculates the completeness score of each pathway and draws the heatmap using the most complete pathways.

NOTE: this sample dataset will not detect all amino acid biosynthesis pathways as the database has been significantly reduced. An example of the output for the reciprocal best hits protocol is available as total_table_biosynthesis_rechit_eval10-40.tsv, while the output for the clustering protocol is available as total_table_biosynthesis_mmseq_c08_eval10-40.tsv in the "examples" directory.


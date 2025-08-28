# MassiveOutsourcing

The repository contains the sample data and the analysis scripts accompanying the manuscript Kasalo, N., Domazet-Lošo, M., Domazet-Lošo, T. **Massive outsourcing of energetically costly amino acids at the origin of animals**.

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
1.

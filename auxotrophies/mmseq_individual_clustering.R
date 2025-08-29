library(readr)
library(plyr)
library(dplyr)
library(stringr)
library(tidyverse)
library(openxlsx)
library(purrr)
library(fuzzyjoin)
library(ggtree)
library(treeio)
library(tidytree)
library(KEGGREST)
library(Biostrings)
library(data.table)
library(fs)

#Continue
#master_table = read_tsv("mmseq_biosynthesis_temp.tsv")
#done = master_table$sseqid
#done = gsub(".*tx|", "", done)
#done = gsub("|", "", done, fixed = TRUE)
#done = paste(done, ".faa", sep = "")
#done = unique(done)

#or From beginning
done = list()

#Define input files
#Fasta file containing the referenca database of enzymes
ref_db = "AA_enzymes_eukarya_sample.faa"
#ref_db = "/Users/nikokasalo/Desktop/Legionella/VF/VFDB_setB_pro.fas"
#Folder containing the individual species proteome fasta files, with headers in the format proteinID|taxID
species_folder = "DB_2024_final"
#species_folder = "mmseq_reduced_further"

#Get a list of files in the folder of individual proteomes
file_list = list.files(species_folder)

#Get the position of the directory in which the script is run
basal_directory = getwd()

#Loop through all files and perform clustering on each individaully, assembling a master table in process
file_num = 1
for (file in file_list) {
  
  if (file %in% done) {
    file_num = file_num + 1
    next
  }
  
  print(file_num)
  #Create a temporary folder for clustering
  dir.create("mmseq_tmp")
  
  #Copy reference database and the current species file to the tmp folder
  file.copy(from = ref_db, to = "mmseq_tmp")
  file.copy(from = paste(species_folder, "/",file, sep = ""), to = "mmseq_tmp")
  
  #Go into the tmp directory
  setwd("mmseq_tmp")
  
  #Concatenate the two fasta files into one named "all.faa"
  system2(command = "cat", args = "* > all.faa")
  
  #Run the mmseqs commands in sequence
  system2(command = "mmseqs", args = "createdb all.faa DB")
  system2(command = "mmseqs", args = "cluster DB DB_clu tmp -c 0.8 --cov-mode 0 --cluster-reassign -e 10E-40")
  system2(command = "mmseqs", args = "createtsv DB DB DB_clu DB_clu.tsv")
  
  #Open the result tsv
  diamond_all = read_tsv("DB_clu.tsv", col_names = FALSE)
  colnames(diamond_all) = c("representative", "member")
  
  #Make two columns, targets: includes only the labels of genes in the target database, genes: query database
  diamond_all$targets = ifelse(grepl("|", diamond_all$member, fixed = TRUE), diamond_all$member, "0")
  diamond_all$genes = ifelse(grepl("|", diamond_all$member, fixed = TRUE), "0", diamond_all$member)
  
  #Make lists of the previous two columns, grouped by representative. This way, we get a list of all queries and targets per each cluster
  diamond_all = diamond_all %>% group_by(representative) %>% mutate(target_list = list(targets))
  diamond_all = diamond_all %>% group_by(representative) %>% mutate(gene_list = list(genes))
  
  #Keep only the last two columns and remove duplicate entries
  diamond_all = diamond_all[5:6]
  diamond_all = unique(diamond_all)
  
  #Unnest the target lists so that each target gets assigned a list of genes it was grouped with
  diamond_all = diamond_all %>% unnest(target_list)
  #Delete the rows that do not contain a target gene
  diamond_all = subset(diamond_all, target_list != "0")
  gc()
  
  #In case of a large table, split into multiple smaller tables
  n = 20000
  df = split(diamond_all, ceiling(seq(nrow(diamond_all))/n))
  length(df)
  
  #Go through subtables and unnest gene lists to get the final table of pairwise connections between genes and the reference database enzymes
  position = 1
  for (i in df) {
    print(position)
    i = i %>% unnest(gene_list)
    i = subset(i, gene_list != "0")
    
    if (position == 1) {
      result_table = i
    } else {
      result_table = rbind(result_table, i)
      print(nrow(result_table))
    }
    position = position + 1
  }
  
  #The final table of pairwise connections
  diamond_all = result_table
  colnames(diamond_all) = c("sseqid", "gene_name")
  gc()
  
  #Add the resulting pairwise connections to a single master table (adds other results from other loop iterations)
  if (file_num == 1) {
    master_table = diamond_all
  } else {
    master_table = rbind(master_table, diamond_all)
    print(nrow(result_table))
  }
  
  #Return to top directory
  setwd(basal_directory)
  
  #Save temporary file
  write_tsv(master_table, "mmseq_biosynthesis_temp.tsv")
  
  #Delete the tmp folder before next iteration
  dir_delete("mmseq_tmp")
  
  #Iterate position
  file_num = file_num + 1
}

#Write the result table for all species together
write_tsv(master_table, "mmseq_biosynthesis_final_08_e-40.tsv")



#Reciprocal best hits####

#Continue
#master_table = read_tsv("mmseq_biosynthesis_reciprocal_temp.tsv")
#done = master_table$qseqid
#done = gsub(".*tx|", "", done)
#done = gsub("|", "", done, fixed = TRUE)
#done = paste(done, ".faa", sep = "")
#done = unique(done)

#or From beginning
done = list()

#Define input files
#Fasta file containing the referenca database of enzymes
ref_db = "AA_enzymes_eukarya_sample.faa"
#ref_db = "/Users/nikokasalo/Desktop/Legionella/VF/VFDB_setB_pro.fas"
#Folder containing the individual species proteome fasta files, with headers in the format proteinID|taxID
species_folder = "DB_2024_final"
#species_folder = "mmseq_reduced_further"

#Get a list of files in the folder of individual proteomes
file_list = list.files(species_folder)


#Get the position of the directory in which the script is run
basal_directory = getwd()

#Loop through all files and perform clustering on each individaully, assembling a master table in process
file_num = 1
for (file in file_list) {
  
  if (file %in% done) {
    file_num = file_num + 1
    next
  }
  
  print(file_num)
  #Create a temporary folder for clustering
  dir.create("mmseq_tmp")
  
  #Copy reference database and the current species file to the tmp folder
  file.copy(from = ref_db, to = "mmseq_tmp")
  file.copy(from = paste(species_folder, "/",file, sep = ""), to = "mmseq_tmp")
  
  #Go into the tmp directory
  setwd("mmseq_tmp")
  
  #Concatenate the two fasta files into one named "all.faa"
  #system2(command = "cat", args = "* > all.faa")
  
  #Run the mmseqs commands in sequence
  #system2(command = "mmseqs", args = "createdb all.faa DB")
  #system2(command = "mmseqs", args = "cluster DB DB_clu tmp -c 0.8 --cov-mode 0 --cluster-reassign")
  #system2(command = "mmseqs", args = "createtsv DB DB DB_clu DB_clu.tsv")
  system2(command = "mmseqs", args = paste("easy-rbh " , file, " ", ref_db, " ABrbh tmp", sep = ""))
  
  #Open the result tsv
  diamond_all = read_tsv("ABrbh", col_names = FALSE)
  #colnames(diamond_all) = c("representative", "member")
  colnames(diamond_all) = c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
  
  
  #Add the resulting pairwise connections to a single master table (adds other results from other loop iterations)
  if (file_num == 1) {
    master_table = diamond_all
  } else {
    master_table = rbind(master_table, diamond_all)
    print(nrow(master_table))
  }
  
  #Return to top directory
  setwd(basal_directory)
  
  #Save temporary file
  write_tsv(master_table, "mmseq_biosynthesis_reciprocal_temp.tsv")
  
  #Delete the tmp folder before next iteration
  dir_delete("mmseq_tmp")
  
  #Iterate position
  file_num = file_num + 1
}

#Write the result table for all species together
write_tsv(master_table, "mmseq_biosynthesis_reciprocal_final.tsv")
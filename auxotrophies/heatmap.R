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

#Using the reciprocal best hit input
diamond_all = read_tsv("mmseq_biosynthesis_reciprocal_final.tsv", col_names = TRUE)
diamond_all = subset(diamond_all, evalue <= 10E-60) 
diamond_all = diamond_all[,c(1:2)]
colnames(diamond_all) = c("sseqid","gene_name")


#Using the mmseqs clustering input
master_table = read_tsv("mmseq_biosynthesis_final_08_e-40.tsv")
diamond_all = master_table
colnames(diamond_all) = c("sseqid", "gene_name")
gc()


#Loading the reference database of enzymes
diamond = read_tsv("pathways_all_KEGG_organisms_eukarya.tsv")
diamond = subset(diamond, select = c(organism, gene_name.2, Reaction.code.x,Pathway,Reaction.code.y,Group_enzyme,Group,Expected))
diamond$gene_name = paste(diamond$organism, ":", diamond$gene_name.2, sep ="")

#Separating genes into focal species

#Add the taxID column
diamond_all$taxID = gsub(".*tx\\|","",diamond_all$sseqid)
diamond_all$taxID = gsub("\\|.*","",diamond_all$taxID)
#diamond_all$taxID = gsub(".*\\|","",diamond_all$sseqid)

taxIDs = unique(diamond_all$taxID)

dir.create("mmseq_all_species_biosynthesis")

pos = 1
for (sp in taxIDs) {
  print(pos)
  sp_table = subset(diamond_all, taxID == sp)
  sp_table = join(sp_table, diamond, by = "gene_name", type = "left")
  sp_table = sp_table[!duplicated(sp_table[c("Reaction.code.x", "Pathway")]),]
  sp_table = subset(sp_table, select = -c(organism, gene_name.2, sseqid))
  write_tsv(sp_table, paste("mmseq_all_species_biosynthesis/", sp, ".enzymes", sep = ""))
  if (pos == 1) {
    total_table = sp_table
  } else {
    total_table = rbind(total_table, sp_table)
  }
  pos = pos + 1
}

#diamond = join(diamond, diamond_all, by = "gene_name", type = "right")

rm(diamond_all)
gc()

total_table = na.omit(total_table)
write_tsv(total_table, "total_table_biosynthesis_rechit_eval10-40.tsv")
total_table = read_tsv("total_table_biosynthesis_rechit_eval10-40.tsv")
diamond = total_table


#Load pathway data
pathways = read.xlsx("Biosynthethic_pathways.xlsx", "All")

#REMOVING ISOLEUCINE_5 TO AVOID CONTAMINATION OF RESULTS - THIS PATH IS ONLY IN ANAEROBIC BACTERIA
pathways = subset(pathways, Group != "Isoleucine_5")

#Add regex symbols for start and end to avoid multiple matches of different enzymes
pathways$Reaction.code = paste("^", pathways$Reaction.code, "$", sep = "")

path_name = pathways$Pathway
enzyme_name = pathways$Reaction.code

#Load the species names-taxID mapping file
#names = read.table("plants_ncbi_dataset.txt", sep = "\t", header = FALSE)
names = read_tsv("names_2023.txt", col_names = FALSE)

colnames(names) = c("taxID", "species")


#diamond_cleaned = subset(diamond, select = c(qseqid, taxID))
#colnames(diamond_cleaned) = c("gene_name", "taxID")

#diamond_cleaned = join(diamond_cleaned, mbr_table[,c(1,3)], by = "gene_name", type = "left")
diamond = join(diamond, names, by = "taxID", type = "left")
#diamond$species = diamond$taxID
diamond$taxID = 1
names(diamond)[names(diamond) == 'taxID'] <- 'presence'

diamond = subset(diamond, select = -c(gene_name))


#write_tsv(diamond, "diamond_checkpoint1.tsv")
#diamond = read_tsv("diamond_checkpoint1.tsv")


diamond = diamond %>% pivot_wider(names_from = Pathway, values_from = presence)

# Remove parentheses from species name so that they correspond to the names in the phylogeny
#diamond$species = gsub("\\s*\\([^\\)]+\\)","",as.character(diamond$species))

diamond = as.data.frame(diamond)
#Change the lists of 1 to "1" to denote presence
#diamond = diamond %>% rowwise %>% mutate(across(7:(ncol(diamond)), ~ toString(.)))
#gc()
#diamond = diamond_test1 %>% mutate(across(7:(ncol(diamond)), ~ case_when(nchar(.)>0 ~ 1)))
diamond = diamond %>% mutate(across(7:(ncol(diamond)), ~ fcase((.)!="NA",1,default = NA)))
gc()
diamond[is.na(diamond)] = 0
diamond = as.data.frame(diamond)


# Unifying the pathway enzymes to produce per-AA results

line3 = diamond %>%
  group_by(species) %>%
  summarise(across(where(is.numeric), .fns = sum))

line3 = line3 %>% mutate(across(3:(ncol(line3)), ~ fcase(.>0,1,default = 0)))
line3 = subset(line3, select = -c(Expected))
line3 = as.data.frame(line3)
rownames(line3) = line3[,1]
line3[,1] = NULL

line3_1 = as.data.frame(t(line3))
line3_1$Pathway = rownames(line3_1)

line3_1 = join(line3_1, pathways[1:4], by = "Pathway", type = "left")

line3_1 = line3_1 %>% mutate(across(1:(ncol(line3_1)-4), as.numeric))


line3_data=line3_1 %>%
  group_by(Group_enzyme) %>%
  summarise(across(where(is.numeric), .fns = sum))

line3_data = line3_data %>% mutate_if(is.numeric, ~1 * (. > 0))

line3_data = join(line3_data, unique(pathways[3:4]), by = "Group_enzyme", type = "left")

line3_data=line3_data %>%
  group_by(Group) %>%
  summarise(across(where(is.numeric), .fns = sum))

line3_data = distinct(join(line3_data, pathways[4:5], by = "Group", type = "left"))
line3_data = na.omit(line3_data)

# Dividing every value in every row by its corresponding expected value of hits. The expected value is the number of
# enzymes expected to be found in a functional pathway, i.e. the amount of steps in the biosynthesis pathway without counting
# multiple possible enzymes per reaction.
line3_data = line3_data %>% rowwise() %>% mutate(across(2:ncol(line3_data), ~ . / Expected))

line3_data = subset(line3_data, select = -c(Expected))

line3_data = as.data.frame(t(line3_data))
colnames(line3_data) = line3_data[1,]
line3_data = line3_data[-1,]
line3_data$species = rownames(line3_data)
line3_data = line3_data %>% relocate(species)



line4 = line3_data %>% mutate(across(where(is.numeric), ~ifelse(. > 0.99, 1, .)))
#line4 = line3_data %>% mutate_if(is.numeric, ~1 * (. > 0.99))
line4 = as.data.frame(line4)

rownames(line4) = line4$species
line4 = line4[2:ncol(line4)]
line4 = line4[1:ncol(line4)]



pathways$AAsynthesis = pathways$Group
pathways$AAsynthesis = gsub("_[1-9]", "", pathways$AAsynthesis)

line5 = t(line4)
line5 = as.data.frame(line5)
line5$Group = rownames(line5)

line5 = join(line5, unique(pathways[4:6]), by = "Group", type = "left")
line5 = subset(line5, select = -c(Expected, Group))

#line5=line5 %>%
#group_by(AAsynthesis) %>%
#summarise(across(where(is.numeric), .fns = sum))
line5=line5 %>%
  mutate_at(c(1:ncol(line5)-1), as.numeric) %>%
  group_by(AAsynthesis) %>%
  summarise(across(where(is.numeric), .fns = max))



#line5 = line5 %>% mutate_if(is.numeric, ~1 * (. > 0.99))

line5 = t(line5)
colnames(line5) = line5[1,]
line5 = line5[-1,]
line5 = as.data.frame(line5)


aa_labels = read_tsv("aa_labels.txt", col_names = T)
energies = read_tsv("listAA20_energy_new_all.txt", col_names = T)
energies = as.data.frame(energies)
energies = join(energies, aa_labels, by = "AA", type = "left")
energies = energies %>% arrange(opportunity.high.respiration)

ordered_AAs = unique(energies$Amino_acid)
line5 = line5[ordered_AAs]
ordered_AAs = intersect(colnames(line5), ordered_AAs)
line5 = line5 %>% mutate_all(as.numeric)

line6 = line5
line6$species = rownames(line6)
write_tsv(line6, "total_table_biosynthesis_mmseq_c08_eval10-40.tsv", col_names = TRUE) 

line6 = as.data.frame(read_tsv("total_table_biosynthesis_mmseq_c08_eval10-40.tsv"))
line6$species = gsub("\\s*\\([^\\)]+\\)","",line6$species)
#line6 = as.data.frame(read_tsv("result_mmseq_auxotrophs.tsv"))
rownames(line6) = line6$species
line6$species = NULL

tree = read.phyloxml("phylo_2023.xml")
#Removing parentheses from the tree (2023 version)
tree_labels = tree@phylo$tip.label
tree_labels = gsub("\\s*\\([^\\)]+\\)","",as.character(tree_labels))
tree@phylo$tip.label = tree_labels

#Subset the tree to Holozoa
tree = tree_subset(tree, node = 642, levels_back = 0)


tree_gg = ggtree(tree, size = 1) + geom_tiplab() + geom_nodelab() + geom_text(aes(label=node), hjust=-.9, vjust = 1)
#tree_gg


tree_gg = ggtree(tree, size = 1.2) + geom_tiplab(size = 6) + geom_nodelab(size = 6, hjust = 1, vjust = -0.4, nudge_x = -0.3) #+ geom_text(aes(label=node), hjust=-.3)
tree_gg = ggtree::rotate(tree_gg, 175)

pdf(file = "heatmap_holozoa_new_c08_eval10-40.pdf",
    width = 35, height = 60)
gheatmap(tree_gg, line6, offset = 16, colnames_angle = 90, colnames_offset_y = 0, font.size = 7, hjust = 1) + 
  scale_fill_gradient2(low = "#F05039", mid = "grey", high = "#3363ff", midpoint = 0.5, guide = "colourbar") +
  guides(fill=guide_legend(title="Completeness\nscore (CS)")) + #scale_fill_discrete(labels = c("Absent", "Present")) + 
  theme(legend.position = c(0.05,0.95), legend.key.size = unit(1.7, "cm"), legend.text = element_text(size=22), legend.title = element_text(size= 24)) + 
  scale_y_discrete(expand = expansion(add = c(5,0))) + scale_x_discrete(expand = expansion(add = c(9,1)))
dev.off()

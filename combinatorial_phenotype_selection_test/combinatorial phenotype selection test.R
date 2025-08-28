library(rcompanion)
library(ggplot2)
library(plyr)
library(dplyr)
library(RcppAlgos)
library(data.table)
library(openxlsx)
library(shadowtext)
library(treeio)
library(ggtree)
library(tidytree)
library(tidyr)
library(readr)
library(factoextra)
library(cocor)
library(ggrepel)
library(gridExtra)
library(rstatix)

#Generating all possible permutations of AA essentiality

permutations_10_10 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(10,10))
permutations_11_9 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(11,9))
permutations_12_8 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(12,8))
permutations_13_7 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(13,7))
permutations_14_6 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(14,6))
permutations_15_5 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(15,5))
permutations_16_4 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(16,4))
permutations_17_3 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(17,3))
permutations_18_2 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(18,2))
permutations_19_1 = permuteGeneral(c("Essential AA", "Non-essential AA"), m = 20, freqs = c(19,1))

permutation_tables = list(permutations_10_10, permutations_11_9, permutations_12_8, permutations_13_7, permutations_14_6,
                          permutations_15_5, permutations_16_4, permutations_17_3, permutations_18_2, permutations_19_1)
auxotroph_count = list(10,11,12,13,14,15,16,17,18,19)

#Loading AA data
aa_labels = read_tsv("aa_labels.txt", col_names = T)
energies = read_tsv("listAA20_energy_new_all.txt", col_names = T)
energies = as.data.frame(energies)
energies = join(energies, aa_labels, by = "AA", type = "left")

#Calculating average values for each AA permutation
for (j in 1:10) {
  
  permutations = permutation_tables[[j]]
  auxotrophies_number = auxotroph_count[[j]]
  
  for (i in 1:nrow(permutations)) {

    if (i%%1000 == 0) {
      print(paste(i, " of ", nrow(permutations), sep = ""))
    }
    iteration = energies
    iteration$shuffled = permutations[i,]
    #iteration$shuffled = sample(iteration$D.extended)
    
    iteration = iteration %>%
      group_by(shuffled) %>% mutate(sum_opportunity_high_resp = sum(opportunity.high.respiration),
                                    mean_opportunity_high_resp = mean(opportunity.high.respiration),
                                    median_opportunity_high_resp = median(opportunity.high.respiration),
                                    sum_cost_high_resp = sum(cost.high.respiration),
                                    mean_cost_high_resp = mean(cost.high.respiration),
                                    median_cost_high_resp = median(cost.high.respiration),
                                    
                                    sum_opportunity_low_resp = sum(opportunity.low.respiration),
                                    mean_opportunity_low_resp = mean(opportunity.low.respiration),
                                    median_opportunity_low_resp = median(opportunity.low.respiration),
                                    sum_cost_low_resp = sum(cost.low.respiration),
                                    mean_cost_low_resp = mean(cost.low.respiration),
                                    median_cost_low_resp = median(cost.low.respiration),
                                    
                                    sum_opportunity_fermentation = sum(opportunity.fermentation),
                                    mean_opportunity_fermentation = mean(opportunity.fermentation),
                                    median_opportunity_fermentation = median(opportunity.fermentation),
                                    sum_cost_fermentation = sum(cost.fermentation),
                                    mean_cost_fermentation = mean(cost.fermentation),
                                    median_cost_fermentation = median(cost.fermentation),
                                    
                              
      )
    
    iteration = iteration %>% distinct(shuffled, sum_opportunity_high_resp, mean_opportunity_high_resp, median_opportunity_high_resp, sum_cost_high_resp, mean_cost_high_resp, median_cost_high_resp,
                                       .keep_all = TRUE)
    iteration = as.data.frame(iteration[10:28])
    rownames(iteration) = iteration[,1]
    iteration = iteration[,-1]
    iteration$n = i
    
    if (i == 1) {
      NE_iteration = iteration["Non-essential AA",]
      E_iteration = iteration["Essential AA",]
    }
    
    if (i > 1) {
      NE_iteration = rbind(NE_iteration, iteration["Non-essential AA",])
      E_iteration = rbind(E_iteration, iteration["Essential AA",])
    }
    
  }
  
  write_tsv(E_iteration, paste("permutations/permutations_",auxotrophies_number,"_elements.tsv", sep = ""))
  write_tsv(NE_iteration, paste("permutations/permutations_",20-auxotrophies_number,"_elements.tsv", sep = ""))
  
}

#Load the AA energy and dispensability data
disp = read.xlsx("disp_cost_table.xlsx")

#Calculate the average cost of the auxotrophic AA set observed in nature
actual_oc_9EAAs = mean(subset(disp$opportunity.high.respiration, disp$D.animals == "Essential AA"))
actual_oc_11EAAs = mean(subset(disp$opportunity.high.respiration, disp$D.extended == "Essential AA"))

#Load the appropriate permutation file and calculate its probability density
n = 11
file_name = paste("permutations/permutations_",n,"_elements.tsv", sep = "")
E_iteration = read_tsv(file_name)

#Calculate the frequency of each average value and use that to calculate the probability of each value
n_elements_oc_high = as.data.frame(table(E_iteration$mean_opportunity_high_resp))
n_elements_oc_high$Var2 = round(as.numeric(levels(n_elements_oc_high$Var1))[n_elements_oc_high$Var1], digits = 2)
n_elements_oc_high$probability = n_elements_oc_high$Freq / sum(n_elements_oc_high$Freq)

#Calculate the p-value of obtaining the identical or higher average value
pPerm = round((sum(n_elements_oc_high[which(n_elements_oc_high$Var2>=round(actual_oc_11EAAs, digits=2)), 4])), digits = 5)


# drawing the graph - PROBABILITY DENSITY FUNCTION - HIGH RESPIRATION - 11 EAAs
#pdf("Figures_rev/11_savings_PDF_high.pdf", width = 9, height = 7)
color_vector = c("black", "black", "#F05039","#F05039")
ggplot(data = n_elements_oc_high, aes(x = Var2, y=probability*100)) + 
  geom_bar(stat="identity", aes(fill = Var2 >= round(actual_oc_11EAAs, digits = 2))) +
  scale_fill_manual(values = c("#6d00c2", "#F05039")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text=element_text(size=20, color = "black"),
        legend.position = "none",
        axis.title = element_text(size=24),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color = color_vector)) +
  annotate("rect", xmin = actual_oc_11EAAs, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.2, fill = "#F05039") +
  #annotate("text", x=actual_oc-2, y=0.17,
  #label= paste(format(actual_oc, format = "e", digits = 4)), size = 8, color = "black") +
  scale_x_continuous(breaks = c(min(n_elements_oc_high$Var2), round(median(n_elements_oc_high$Var2), digits = 2), 
                                round(actual_oc_11EAAs, digits = 2), max(n_elements_oc_high$Var2))) + 
  geom_vline(aes(xintercept = actual_oc_11EAAs), color = "#F05039", linetype = "dashed", size = 1) +
  annotate("text", x=actual_oc_11EAAs, y=max(n_elements_oc_high$probability)*100, hjust = 1,
           label= paste("p =", format(sum(n_elements_oc_high[which(n_elements_oc_high$Var2>=round(actual_oc_11EAAs, digits=2)), 4]), 
                                      format = "e", digits = 2)), size = 8, color = "black") +
  annotate("text", x = min(n_elements_oc_high$Var2), y = max(n_elements_oc_high$probability)*100,
           label = paste("high respiration (",n," auxotrophic AAs)", sep = ""), size = 6.5, color = "black", hjust = 0) +
  xlab(expression(bold("Average opportunity cost (~P)"))) +
  ylab(expression(bold("Proportion of permutations (%)"))) +
  #ggtitle("All permutations (11 EAAs)") +
  theme(plot.title = element_text(size = 20)) + 
  scale_y_continuous(expand = expansion(mult = c(0, .1)))
  #ggtitle(paste(i, " (",unique(plot_table$wide_group),")", sep =""))
#dev.off()
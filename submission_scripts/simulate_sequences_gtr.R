
all_molecular_data_size <- 100000
medium_molecular_data_size <- 10000
small_molecular_data_size <- 1000

##########GENERATE_INDIVIDUAL_HETEROGENEITY/MOLECULAR_EVOLUTIONARY_RATE_SEQUENCES_WITH_DIFFEREN_QUANTITIES_OF_DATA

all_sequences_all_molecular_data <- simSeq(outgroup_clock_tree, l = all_molecular_data_size, type = "DNA", rate = 1, Q=1:6)
all_sequences_medium_molecular_data <- subset(all_sequences_all_molecular_data, select = seq(1, medium_molecular_data_size, 1), site.pattern = FALSE)
all_sequences_small_molecular_data <- subset(all_sequences_all_molecular_data, select = seq(1, small_molecular_data_size, 1), site.pattern = FALSE)

###########GENERATE_DIFFERENT_INDIVIDUAL_HETEROGENEITY/MOLECULAR_EVOLUTIONARY_RATE_DATASETS

three_taxon_names <- c(paste(seq(1, length(three_taxon_trees), 1), ".nexus", sep = ""))
three_taxon_tree_names <- c(paste(seq(1, length(three_taxon_trees), 1), ".tre", sep = ""))

setwd(directories[1])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_all_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

setwd("C:/Users/some3165/Desktop/RATE_TIME_FINAL_2/PHYLOGENETIC_INFERENCE")
setwd(directories[2])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_medium_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

setwd("C:/Users/some3165/Desktop/RATE_TIME_FINAL_2/PHYLOGENETIC_INFERENCE")
setwd(directories[3])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_small_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

######

three_hundered_taxon_tree_names <- "1.tre"
three_hundered_taxon_entire_tree_names <- "1_three_hundered_taxon_entire_tree.tre"
three_hundered_taxon_relaxed_clock_tree_names <- "1_rate.tre"

setwd("C:/Users/some3165/Desktop/RATE_TIME_FINAL_2/PHYLOGENETIC_INFERENCE")
write.tree(entire_tree, three_hundered_taxon_entire_tree_names)
write.tree(outgroup_clock_tree, three_hundered_taxon_relaxed_clock_tree_names)

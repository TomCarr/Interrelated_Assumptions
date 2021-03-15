
all_molecular_data_size <- 100000
medium_molecular_data_size <- 10000
small_molecular_data_size <- 1000

##########GENERATE_CLOCK_SEQUENCES_WITH_DIFFEREN_QUANTITIES_OF_DATA

all_sequences_all_molecular_data <- simSeq(clock_tree, l = all_molecular_data_size, type = "DNA", rate = 1)
all_sequences_medium_molecular_data <- subset(all_sequences_all_molecular_data, select = seq(1, medium_molecular_data_size, 1), site.pattern = FALSE)
all_sequences_small_molecular_data <- subset(all_sequences_all_molecular_data, select = seq(1, small_molecular_data_size, 1), site.pattern = FALSE)

###########GENERATE_DIFFERENT_INDIVIDUAL_HETEROGENEITY_DATASETS

three_taxon_names <- c(paste(seq(1, length(three_taxon_trees), 1), ".nexus", sep = ""))
three_taxon_tree_names <- c(paste(seq(1, length(three_taxon_trees), 1), ".tre", sep = ""))

setwd(directories[1])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_all_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")
setwd(directories[2])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_medium_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")
setwd(directories[3])
for (i in 1:length(three_taxon_trees)){
write.phyDat(subset(all_sequences_small_molecular_data, three_taxon_trees[[i]]$tip.label), file = three_taxon_names[[i]], format = "nexus")
write.tree(three_taxon_trees[[i]], three_taxon_tree_names[[i]])
}

######

thirty_taxon_small_names <- c(paste(seq(1, length(thirty_taxon_trees), 1), ".nexus", sep = "")) 
thirty_taxon_tree_names <- c(paste(seq(1, length(thirty_taxon_trees), 1), ".tre", sep = ""))
three_hundered_taxon_small_names <- "1.nexus"
three_hundered_taxon_tree_names <- "1.tre"
three_hundered_taxon_entire_tree_names <- "1_three_hundered_taxon_entire_tree.tre"
three_hundered_taxon_clock_tree_names <- "1_clock.tre"

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")
setwd(directories[4])
for (i in 1:length(thirty_taxon_trees)){
write.phyDat(subset(all_sequences_small_molecular_data, thirty_taxon_trees[[i]]$tip.label), file = thirty_taxon_small_names[[i]], format = "nexus")
write.tree(thirty_taxon_trees[[i]], thirty_taxon_tree_names[[i]])
}

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")
setwd(directories[5])
write.phyDat(subset(all_sequences_small_molecular_data, extant_tree$tip.label), file = three_hundered_taxon_small_names, format = "nexus")
write.tree(extant_tree, three_hundered_taxon_tree_names)
write.tree(entire_tree, three_hundered_taxon_entire_tree_names)
write.tree(clock_tree, three_hundered_taxon_clock_tree_names)

######

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")
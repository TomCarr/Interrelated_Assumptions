library(phangorn)
library(devtools)
library(NELSI)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(ggplot2)
library(Rmisc)


######EXTRACT_FOSSIL_CALIBRATIONS_FOR_300_TAXON_TREE########
############################################################

fossils_ages_file <- "1_r_output_calibration_fossil_ages.txt"
oldest_fossils_amongst_all_fossils_ages_file <- "1_r_output_oldest_fossil_ages.txt"
taxon_sets_file <- "1_taxon_sets.txt"
calibrated_file <- "1_contains_calibrations.txt"

######EXTRACT_FOSSIL_CALIBRATIONS#####

node_vector <- length(extant_tree$tip.label) + seq(1, length(extant_tree$tip.label)-1, 1)
tree_nested_clades <- vector("list", length(node_vector))

for (i in 1:length(tree_nested_clades)){
tree_nested_clades[[i]] <-  extract.clade(tree_with_fossil_tips, findMRCA(tree_with_fossil_tips, extract.clade(extant_tree, node_vector[[i]])$tip.label, type = "node")) # find MRCAs of extant taxa in fossil tree
}

tree_nested_clades_fossils <- vector("list", length(node_vector))
for (j in 1:length(tree_nested_clades)){
for (i in 1:length(tree_nested_clades[[j]]$tip.label)){
if (length(unlist(strsplit(tree_nested_clades[[j]]$tip.label[[i]], split = "t"))) == 1){
tree_nested_clades_fossils[[j]] <- append(tree_nested_clades_fossils[[j]], tree_nested_clades[[j]]$tip.label[[i]])
}
}
}

######GET_AGES_OF_EXTRACTED_FOSSILS#####

tree_nested_clades_fossils_ages <- vector("list", length(node_vector))

for (j in 1:length(tree_nested_clades_fossils_ages)){
if (length(tree_nested_clades_fossils[[j]]) > 0){
for (i in 1:length(tree_nested_clades_fossils[[j]])){
tree_nested_clades_fossils_ages[[j]] <- append(tree_nested_clades_fossils_ages[[j]], unlist(absolute_fossil_ages)[[sapply(tree_nested_clades_fossils[[j]][[i]], as.numeric)[[1]]]])
}
}
}

######GET_OLDEST_FOSSIL_FOR_EACH_CLADE#####

oldest_fossils <- vector("list", length(node_vector))

for (i in 1:length(tree_nested_clades_fossils_ages)){
if (length(tree_nested_clades_fossils_ages[[i]]) > 0 ){
oldest_fossils[[i]] <- max(tree_nested_clades_fossils_ages[[i]])
}
}

#######REMOVE_FOSSILS_WITH_OLDER_FOSSILS_IN_DESCENDENT_CLADES#####

older_descendant_node_calibration <- vector(mode="numeric", length=0)
descendants_of_calibrated_nodes <- vector("list", length(node_vector))

for (i in 1:length(oldest_fossils)){
descendants_of_calibrated_nodes[[i]] <- getDescendants(extant_tree, length(extant_tree$tip.label) + i)
}

tip_label_removal_vector <- vector("list", length(node_vector))
for (j in 1:length(descendants_of_calibrated_nodes)){
for (i in 1:length(descendants_of_calibrated_nodes[[j]])){
if (descendants_of_calibrated_nodes[[j]][[i]] <= length(extant_tree$tip.label)){
tip_label_removal_vector[[j]] <- append(tip_label_removal_vector[[j]], i)
}
}
}

for (i in 1:length(descendants_of_calibrated_nodes)){
descendants_of_calibrated_nodes[[i]] <- descendants_of_calibrated_nodes[[i]][-tip_label_removal_vector[[i]]]
}

for (j in 1:length(descendants_of_calibrated_nodes)){
if (length(descendants_of_calibrated_nodes[[j]]) > 0 ){
for (i in 1:length(descendants_of_calibrated_nodes[[j]])){
descendants_of_calibrated_nodes[[j]][[i]] <- descendants_of_calibrated_nodes[[j]][[i]] - (length(extant_tree$tip.label))
}
}
}

older_descendant <- vector("list", length(node_vector))
for (j in 1:length(descendants_of_calibrated_nodes)){
if (length(descendants_of_calibrated_nodes[[j]]) > 0){
for (i in 1:length(descendants_of_calibrated_nodes[[j]])){
if (length(oldest_fossils[[descendants_of_calibrated_nodes[[j]][[i]]]]) > 0){
if (oldest_fossils[[descendants_of_calibrated_nodes[[j]][[i]]]] > oldest_fossils[[j]]){
older_descendant[[j]] <- append(older_descendant[[j]], 1)
}
}
}
}
}

for (i in 1:length(oldest_fossils)){
if (length(older_descendant[[i]]) < 0){
oldest_fossils[[i]] <- vector(mode="numeric", length=0)
}
}

######GET_TAXON_SETS_FOR_CALIBRATED_CLADES#####

calibrated_clade_taxon_sets <- vector("list", length(oldest_fossils))
for (i in 1:length(oldest_fossils)){
if (length(oldest_fossils[[i]]) > 0){
calibrated_clade_taxon_sets[[i]] <- append(calibrated_clade_taxon_sets[[i]], extract.clade(extant_tree, node_vector[[i]])$tip.label)
}
}

######REMOVE_RECORDS_OF_CLADES_WITH_NO_FOSSIL_CALIBRATIONS#####

zero_length_node_calibration_removal <- vector(mode="numeric", length=0)
for (i in 1:length(oldest_fossils)){
if (length(oldest_fossils[[i]]) == 0){
zero_length_node_calibration_removal <- append(zero_length_node_calibration_removal, i)
}
}

if (length(zero_length_node_calibration_removal) > 0){
oldest_fossils <- oldest_fossils[-zero_length_node_calibration_removal]
calibrated_clade_taxon_sets <- calibrated_clade_taxon_sets[-zero_length_node_calibration_removal]
}

######BUT_CONDITION_SO_VECTOR_ACTUALLY_EXISTS

if (length(calibrated_clade_taxon_sets) == 0){
oldest_fossils[[length(oldest_fossils) + 1]] <- 0
calibrated_clade_taxon_sets[[length(calibrated_clade_taxon_sets) + 1]] <- c("t", "t", "t")
}

######GET_OLDEST_AMONGST_ALL_FOSSILS

oldest_amongst_all_fossils <- vector(mode="numeric", length=0)
for (i in 1:length(unlist(oldest_fossils))){
oldest_amongst_all_fossils <- max(unlist(oldest_fossils))
}

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE/three_hundered_taxon_small/1")

######WRITE_FOSSILS_TO_FILE

if (oldest_fossils[[1]] == 0){
sink(fossils_ages_file)
print(oldest_fossils)
sink()
sink(taxon_sets_file)
print(calibrated_clade_taxon_sets)
sink()
sink(oldest_fossils_amongst_all_fossils_ages_file)
print(oldest_amongst_all_fossils)
sink()
calibrate <- 0
sink(calibrated_file)
print(calibrate)
sink()
} else {
sink(fossils_ages_file)
print(oldest_fossils)
sink()
sink(taxon_sets_file)
print(calibrated_clade_taxon_sets)
sink()
sink(oldest_fossils_amongst_all_fossils_ages_file)
print(oldest_amongst_all_fossils)
sink()
calibrate <- 1
sink(calibrated_file)
print(calibrate)
sink()
}

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/DIVERSIFICATION_RATE/ENTIRE_TREE")

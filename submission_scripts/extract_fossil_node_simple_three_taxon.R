library(phangorn)
library(devtools)
library(NELSI)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(ggplot2)
library(Rmisc)

######EXTRACT_FOSSIL_CALIBRATIONS_FOR_THREE_TAXON_TREES#####
############################################################

fossils_ages_files <- paste(seq(1, length(three_taxon_trees), 1), "_r_output_calibration_fossil_ages.txt", sep = "")
oldest_fossils_amongst_all_fossils_ages_files <- paste(seq(1, length(three_taxon_trees), 1), "_r_output_oldest_fossil_ages.txt", sep = "")
taxon_sets_files <- paste(seq(1, length(three_taxon_trees), 1), "_taxon_sets.txt", sep = "")
calibrated_files <- paste(seq(1, length(three_taxon_trees), 1), "_contains_calibrations.txt", sep = "")

for (a in 1:length(three_taxon_trees)){

three_taxon_tree <- three_taxon_trees[[a]]

######EXTRACT_FOSSIL_CALIBRATIONS#####

node_vector <- length(three_taxon_tree$tip.label) + seq(1, length(three_taxon_tree$tip.label)-1, 1)
three_taxon_tree_nested_clades <- vector("list", length(node_vector))

for (i in 1:length(three_taxon_tree_nested_clades)){
three_taxon_tree_nested_clades[[i]] <- extract.clade(tree_with_fossil_tips, findMRCA(tree_with_fossil_tips, extract.clade(three_taxon_tree, node_vector[[i]])$tip.label, type = "node"))
}

three_taxon_tree_nested_clades_fossils <- vector("list", length(node_vector))
for (j in 1:length(three_taxon_tree_nested_clades)){
for (i in 1:length(three_taxon_tree_nested_clades[[j]]$tip.label)){
if (length(unlist(strsplit(three_taxon_tree_nested_clades[[j]]$tip.label[[i]], split = "t"))) == 1){
three_taxon_tree_nested_clades_fossils[[j]] <- append(three_taxon_tree_nested_clades_fossils[[j]], three_taxon_tree_nested_clades[[j]]$tip.label[[i]])
}
}
}

######GET_AGES_OF_EXTRACTED_FOSSILS#####

three_taxon_tree_nested_clades_fossils_ages <- vector("list", length(node_vector))

for (j in 1:length(three_taxon_tree_nested_clades_fossils_ages)){
if (length(three_taxon_tree_nested_clades_fossils[[j]]) > 0){
for (i in 1:length(three_taxon_tree_nested_clades_fossils[[j]])){
three_taxon_tree_nested_clades_fossils_ages[[j]] <- append(three_taxon_tree_nested_clades_fossils_ages[[j]], unlist(absolute_fossil_ages)[[sapply(three_taxon_tree_nested_clades_fossils[[j]][[i]], as.numeric)[[1]]]])
}
}
}

######GET_OLDEST_FOSSIL_FOR_EACH_CLADE#####

oldest_fossils <- vector("list", length(node_vector))

for (i in 1:length(three_taxon_tree_nested_clades_fossils_ages)){
if (length(three_taxon_tree_nested_clades_fossils_ages[[i]]) > 0 ){
oldest_fossils[[i]] <- max(three_taxon_tree_nested_clades_fossils_ages[[i]])
}
}

#######REMOVE_FOSSILS_WITH_OLDER_FOSSILS_IN_DESCENDENT_CLADES#####

older_descendant_node_calibration <- vector(mode="numeric", length=0)
descendants_of_calibrated_nodes <- vector("list", length(node_vector))

for (i in 1:length(oldest_fossils)){
descendants_of_calibrated_nodes[[i]] <- getDescendants(three_taxon_tree, length(three_taxon_tree$tip.label) + i)
}

tip_label_removal_vector <- vector("list", length(node_vector))
for (j in 1:length(descendants_of_calibrated_nodes)){
for (i in 1:length(descendants_of_calibrated_nodes[[j]])){
if (descendants_of_calibrated_nodes[[j]][[i]] <= length(three_taxon_tree$tip.label)){
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
descendants_of_calibrated_nodes[[j]][[i]] <- descendants_of_calibrated_nodes[[j]][[i]] - (length(three_taxon_tree$tip.label))
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
calibrated_clade_taxon_sets[[i]] <- append(calibrated_clade_taxon_sets[[i]], extract.clade(three_taxon_tree, node_vector[[i]])$tip.label)
}
}

######REMOVE_RECORDS_OF_CLADES_WITH_NO_FOSSIL_CALIBRATIONS

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

setwd("C:/Users/some3165/Documents/RATE_TIME_PROBLEM_FINAL/INDIVIDUAL_HETEROGENEITY/MOLECULAR_EVOLUTIONARY_RATE_WITH_YULE")
setwd(directories[1])

######WRITE_FOSSILS_TO_FILE

if (oldest_fossils[[1]] == 0){
sink(fossils_ages_files[[a]])
print(oldest_fossils)
sink()
sink(taxon_sets_files[[a]])
print(calibrated_clade_taxon_sets)
sink()
sink(oldest_fossils_amongst_all_fossils_ages_files[[a]])
print(oldest_amongst_all_fossils)
sink()
calibrate <- 0
sink(calibrated_files[[a]])
print(calibrate)
sink()
} else {
sink(fossils_ages_files[[a]])
print(oldest_fossils)
sink()
sink(taxon_sets_files[[a]])
print(calibrated_clade_taxon_sets)
sink()
sink(oldest_fossils_amongst_all_fossils_ages_files[[a]])
print(oldest_amongst_all_fossils)
sink()
calibrate <- 1
sink(calibrated_files[[a]])
print(calibrate)
sink()
}
}


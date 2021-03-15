library(phangorn)
library(devtools)
library(NELSI)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(ggplot2)
library(Rmisc)

rate <- 5

entire_tree_branching_times <- round(branching.times(entire_tree), digits = 5) + (max(node.depth.edgelength(entire_tree)) - round(branching.times(entire_tree)[[1]], digits = 5))

######CONSTRUCT_FOSSIL_SIMULATION_DATA_FRAME#####
#################################################

######CONSTRUCT_ORDERED_FOSSIL_NODES#####

fossil_simulation_data_frame <- data.frame(entire_tree$edge[,1], entire_tree$edge[,2], entire_tree$edge.length)
fossil_simulation_data_frame <- fossil_simulation_data_frame[order(fossil_simulation_data_frame[,1]),]


######GET_BRANCHING_TIMES_OF_ANCESTRAL_NODES#####

ancestral_branching_times <- vector("list", length(entire_tree_branching_times))
for (i in 1:length(entire_tree_branching_times)){
ancestral_branching_times[[i]] <- c(rep(entire_tree_branching_times[i], 2))
}
ancestral_branching_times <- unlist(ancestral_branching_times)

fossil_simulation_data_frame <- cbind(fossil_simulation_data_frame, ancestral_branching_times)


######SIMULATE_FOSSIL_OCCURENCES_ACCORDING_TO_A_GIVEN_RATE_FOR_EACH_BRANCH_IN_THE_PHYLOGENY#####

fossil_string <- vector("list", nrow(fossil_simulation_data_frame))
for (i in 1:length(fossil_string)){
fossil_string[[i]] <- rexp(1, rate)
}

for (i in 1:length(fossil_string)){
while (fossil_string[[i]][[length(fossil_string[[i]])]] < fossil_simulation_data_frame[,3][[i]]){
fossil_string[[i]] <- append(fossil_string[[i]], fossil_string[[i]][[length(fossil_string[[i]])]] + rexp(1, rate))
}
}

for (j in 1:length(fossil_string)){
for (i in 1:length(fossil_string[[j]])){
if (fossil_string[[j]][[i]] > fossil_simulation_data_frame[,3][[j]]){
fossil_string[[j]][[i]] <- 1000000
} else {
fossil_string[[j]][[i]] <- fossil_string[[j]][[i]]
}
}
}

######REMOVE_FOSSILS_SIMULATED_AFTER_THE_END_OF_A_BRANCH#####

removal_fossils <- vector("list", length(fossil_string))
for (j in 1:length(fossil_string)){
for (i in 1:length(fossil_string[[j]])){
if (fossil_string[[j]][[i]] == 1000000){
removal_fossils[[j]] <- append(removal_fossils[[j]], i)
}
}
}

for(j in 1:length(fossil_string)){
if (length(removal_fossils[[j]]) > 0){
fossil_string[[j]] <- fossil_string[[j]][-removal_fossils[[j]]]
}
}

######GET_ABSOLUTE_FOSSIL_AGES#####

absolute_fossil_ages <- vector("list", nrow(fossil_simulation_data_frame))
for (j in 1:length(absolute_fossil_ages)){
if (length(fossil_string[[j]]) > 0){
for (i in 1:length(fossil_string[[j]])){
absolute_fossil_ages[[j]][[i]] <- fossil_simulation_data_frame[,4][[j]] - fossil_string[[j]][[i]]
}
}
}

######GET_FOSSIL_BRANCH_LENGTH#####

repeated_fossil_branch_lengths <- vector("list", length(absolute_fossil_ages))
for (i in 1:length(absolute_fossil_ages)){
if (length(absolute_fossil_ages[[i]]) > 0){
repeated_fossil_branch_lengths[[i]] <- rep(fossil_simulation_data_frame[,3][[i]], length(absolute_fossil_ages[[i]]))
}
}

repeated_fossil_branch_lengths <- unlist(repeated_fossil_branch_lengths)

######ADD_SIMULATED_FOSSILS_TO_SIMULATED_TREE#####
##################################################

if (length(unlist(absolute_fossil_ages)) > 0){

######GET_POSITIONS_OF_FOSSILS#####

altered_fossil_nodes_descendant <- vector("list", nrow(fossil_simulation_data_frame))
for (j in 1:nrow(fossil_simulation_data_frame)){
if ((length(absolute_fossil_ages[[j]]) > 0) & (fossil_simulation_data_frame[,2][[j]] > length(entire_tree$tip.label))){
for (i in 1:length(absolute_fossil_ages[[j]])){
altered_fossil_nodes_descendant[[j]][[i]] <- fossil_simulation_data_frame[,2][[j]]
}
} else if (length(absolute_fossil_ages[[j]] > 0) & (fossil_simulation_data_frame[,2][[j]] <= length(entire_tree$tip.label))){
for (i in 1:length(absolute_fossil_ages[[j]])){
altered_fossil_nodes_descendant[[j]][[i]] <- entire_tree$tip.label[fossil_simulation_data_frame[,2][[j]]]
}
}
}

######SIMULATE_FOSSIL_NAMES_AND_BRANCHES#####

fossil_names <- seq(1, length(unlist(absolute_fossil_ages)), 1)
fossil_branches <- vector("list", length(unlist(absolute_fossil_ages)))
for (i in 1:length(fossil_branches)){
fossil_branches[[i]] <- list(edge=matrix(c(2,1),1,2),tip.label=fossil_names[[i]], edge.length=0, Nnode=1)
class(fossil_branches[[i]]) <- "phylo"
}

######SPLIT_FOSSILS_DEPENDING_ON_WHETHER_THEY_ARE_BEING_ATTACHED_TO_TIP_BRANCHES#####

altered_fossil_nodes_descendant_numeric <- vector("list", nrow(fossil_simulation_data_frame))
altered_fossil_nodes_descendant_character <- vector("list", nrow(fossil_simulation_data_frame))
for (i in 1:length(altered_fossil_nodes_descendant)){
if (length(altered_fossil_nodes_descendant[[i]] > 0)){
if (is.numeric(altered_fossil_nodes_descendant[[i]][[1]]) == TRUE){
altered_fossil_nodes_descendant_numeric[[i]] <- altered_fossil_nodes_descendant[[i]]
} else {
altered_fossil_nodes_descendant_character[[i]] <- altered_fossil_nodes_descendant[[i]]
}
}
}

######CONDITION_THE_TWO_ALTERNATIVE_VECTORS#####

for (i in 1:length(altered_fossil_nodes_descendant_numeric)){
if (length(altered_fossil_nodes_descendant_character[[i]]) > 0){
altered_fossil_nodes_descendant_numeric[[i]] <- append(altered_fossil_nodes_descendant_numeric[[i]], rep(1000000, length(altered_fossil_nodes_descendant_character[[i]])))
} else if (length(altered_fossil_nodes_descendant_numeric[[i]]) > 0){
altered_fossil_nodes_descendant_character[[i]] <- append(altered_fossil_nodes_descendant_character[[i]], rep("GAP", length(altered_fossil_nodes_descendant_numeric[[i]])))
}
}

altered_fossil_nodes_descendant_numeric <- unlist(altered_fossil_nodes_descendant_numeric)
altered_fossil_nodes_descendant_character <- unlist(altered_fossil_nodes_descendant_character)

######BIND_FOSSIL_BRANCHES_TO_TREE#####

tree_with_fossil_tips <- entire_tree

for (i in 1:length(fossil_branches)){
if (altered_fossil_nodes_descendant_numeric[[i]] != 1000000){
tree_with_fossil_tips <- bind.tree(tree_with_fossil_tips, 
fossil_branches[[i]],
where = findMRCA(tree_with_fossil_tips, extract.clade(entire_tree, altered_fossil_nodes_descendant_numeric[[i]])$tip.label, type=c("node")),
position = repeated_fossil_branch_lengths[[i]] - unlist(fossil_string)[[i]])
} else {
tree_with_fossil_tips <- bind.tree(tree_with_fossil_tips,
fossil_branches[[i]],
where = which(tree_with_fossil_tips$tip.label==altered_fossil_nodes_descendant_character[[i]]), position = repeated_fossil_branch_lengths[[i]] - unlist(fossil_string)[[i]])
}
}

}

######GENERATE_FOSSIL_TREE_WITHOUT_EXTINCT_TAXA#####

fossil_tree_without_extinct_taxa <- drop.tip(tree_with_fossil_tips, tip = getExtinct(entire_tree))
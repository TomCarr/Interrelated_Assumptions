#' compiles node age estimates from three taxon trees into single list in order of nodes in species tree
#' @export

CompileAgeEstimates <- function(n_three_taxon_trees, species_tree, outputs_dir){

####################################################################
####################################################################
###GET_SPECIES_TREE_DATA############################################
####################################################################
####################################################################

###GET_SPECIES_TREE_CLADES_AND_TWO_SUBSET_HALVES####################

species_tree_clades <- vector("list", length(species_tree$tip.label)-1)
species_tree_clades_subsets <- vector("list", length(species_tree$tip.label)-1) 
for (i in 1:length(species_tree_clades_subsets)){
species_tree_clades_subsets[[i]] <- vector("list", 2)
}

for (i in 1:length(species_tree_clades)){
species_tree_clades[[i]] <- extract.clade(species_tree, i + length(species_tree$tip.label))

for (a in which(species_tree[[1]][,1] == length(species_tree$tip.label) + i)[[1]]){
if (species_tree[[1]][,2][[a]] > length(species_tree$tip.label)){
species_tree_clades_subsets[[i]][[1]] <- extract.clade(species_tree, species_tree[[1]][,2][[a]])$tip.label
} else {
species_tree_clades_subsets[[i]][[1]] <- species_tree$tip.label[[species_tree[[1]][,2][[a]]]]
}
}

for (a in which(species_tree[[1]][,1] == length(species_tree$tip.label) + i)[[2]]){
if (species_tree[[1]][,2][[a]] > length(species_tree$tip.label)){
species_tree_clades_subsets[[i]][[2]] <- extract.clade(species_tree, species_tree[[1]][,2][[a]])$tip.label
} else {
species_tree_clades_subsets[[i]][[1]] <- species_tree$tip.label[[species_tree[[1]][,2][[a]]]]
}
}

}

###################################################################
###################################################################
###READ_IN_THREE_TAXON_TREES_AND_GET_INFERRED_TIPS#################
###################################################################
###################################################################

three_taxon_trees <- vector("list", n_three_taxon_trees)
inferred_tips <- vector("list", n_three_taxon_trees)
for (i in 1:n_three_taxon_trees){
three_taxon_trees[[i]] <- read.tree(paste(paste(outputs_dir, paste("three_taxon_", i, sep=""), sep=""),"/three_taxon.tre", sep=""))
inferred_tips[[i]] <- append(inferred_tips[[i]], extract.clade(three_taxon_trees[[i]], 5)$tip.label)
}

###################################################################
###################################################################
###APPEND_TO_AGE_ESTIMATE_LIST_FOR_SPECIES_TREE####################
###################################################################
###################################################################

clade_age_estimates <<- vector("list", length(species_tree$tip.label)-1)
for (i in 1:length(species_tree_clades_subsets)){
for (a in 1:length(inferred_tips)){
if (((inferred_tips[[a]][[1]] %in% species_tree_clades_subsets[[i]][[1]]) & (inferred_tips[[a]][[2]] %in% species_tree_clades_subsets[[i]][[2]])) 
|| ((inferred_tips[[a]][[2]] %in% species_tree_clades_subsets[[i]][[1]]) & (inferred_tips[[a]][[1]] %in% species_tree_clades_subsets[[i]][[2]]))){
clade_age_estimates[[i]] <<- append(clade_age_estimates[[i]], max(node.depth.edgelength(three_taxon_trees[[a]])) - findMRCA(three_taxon_trees[[a]], inferred_tips[[a]], "height"))
}
}
}

}

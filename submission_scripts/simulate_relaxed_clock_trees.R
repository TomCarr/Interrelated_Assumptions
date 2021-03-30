library(phangorn)
library(devtools)
library(NELSI)
library(TreeSim)
library(phytools)
library(EnvStats)
library(phylobase)
library(ggplot2)
library(Rmisc)

n_tips <- 301
lambda <- 3
mu <- 0
n_reps <- 1

##########

tip_vector <- seq(1, n_tips, 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

##########

entire_tree <- sim.bd.taxa(n_tips, n_reps, lambda, mu, frac = 1, complete = TRUE, stochsampling = FALSE)
entire_tree <- entire_tree[[1]]

clock_tree <- entire_tree
for (i in 1:length(clock_tree$edge.length)){
clock_tree$edge.length[[i]] <- clock_tree$edge.length[[i]]*rlnorm(1, log(0.05), 0.5)
}

tip <- list(edge=matrix(c(2,1),1,2),tip.label="t302", edge.length=0.1,Nnode=1)
class(tip) <- "phylo"
outgroup_clock_tree <- bind.tree(clock_tree,tip,where=302, position = 0.01)

##########GENERATE_THREE_TAXON_SAMPLES

species <- vector("list", length(node_numbers)/2)
sample_nodes <- vector("list", length(node_numbers)/2)
remaining_node_numbers <- node_numbers

for (i in 1:length(sample_nodes)){
if (length(remaining_node_numbers) > 1 ){
sample_nodes[[i]] <- sample(remaining_node_numbers, 2)
remaining_node_numbers <- remaining_node_numbers[! remaining_node_numbers %in% sample_nodes[[i]]]
} else {
sample_nodes[[i]] <- remaining_node_numbers
}
}

write(unlist(sample_nodes), "two_nodes.txt")

for (j in 1:length(sample_nodes)){
for (i in 1:length(sample_nodes[[j]])){
full_extraction <- extract.clade(entire_tree, sample_nodes[[j]][[i]])
if (length(full_extraction$tip.label) > 2){
not_root <- 1
while (not_root == 1){
two_taxon <- drop.tip(full_extraction, full_extraction$tip.label[-sample(seq(1, length(full_extraction$tip.label), 1), 2)])
if (round(max(node.depth.edgelength(two_taxon)), digits = 4) == round(max(node.depth.edgelength(full_extraction)), digits = 4)){
species[[j]] <- append(species[[j]], two_taxon$tip.label)
not_root <- 0
}
}
} else {
two_taxon <- full_extraction
species[[j]] <- append(species[[j]], two_taxon$tip.label)
}
}
print(j)
}

sample_two_nodes <- vector("list", length(species))
for (i in 1:length(species)){
sample_two_nodes[[i]] <- unique(species[[i]])
}

##########THREE_TAXON_TREES

two_node_index <- vector("list", length(sample_two_nodes))
for (j in 1:length(sample_two_nodes)){
for (i in 1:length(sample_two_nodes[[j]])){
two_node_index[[j]] <- append(two_node_index[[j]], which(outgroup_clock_tree$tip.label == sample_two_nodes[[j]][[i]])) 
}
}

for (i in 1:length(sample_two_nodes)){
two_node_index[[i]] <- append(two_node_index[[i]], which(outgroup_clock_tree$tip.label == sample(outgroup_clock_tree$tip.label[-which(outgroup_clock_tree$tip.label %in% extract.clade(clock_tree, findMRCA(clock_tree, two_node_index[[i]], type = "node"))$tip.label)], 1)))
}

three_taxon_trees <- vector("list", length(two_node_index))
for (i in 1:length(three_taxon_trees)){
three_taxon_trees[[i]] <- drop.tip(outgroup_clock_tree, outgroup_clock_tree$tip.label[-two_node_index[[i]]])
}









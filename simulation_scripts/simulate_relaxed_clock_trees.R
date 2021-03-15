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

relaxed_clock_tree <- entire_tree
for (i in 1:length(relaxed_clock_tree$edge.length)){
relaxed_clock_tree$edge.length[[i]] <- relaxed_clock_tree$edge.length[[i]]*rlnorm(1, log(0.05) - ((0.587405^2)/2), sdlog=0.587405)
}

extant_tree <- drop.tip(entire_tree, getExtinct(entire_tree))
root_age <- max(node.depth.edgelength(extant_tree))
write(root_age, "root_age.txt")

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
full_extraction <- extract.clade(extant_tree, sample_nodes[[j]][[i]])
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

##########GENERATE_THIRTY_TAXON_SAMPLES

species <- vector("list", length(node_numbers)/30)
sample_nodes <- vector("list", length(node_numbers)/30)
remaining_node_numbers <- node_numbers

for (i in 1:length(species)){
if (length(remaining_node_numbers) > 0){
sample_nodes[[i]] <- sample(remaining_node_numbers, 30)
remaining_node_numbers <- remaining_node_numbers[! remaining_node_numbers %in% sample_nodes[[i]]]
} else {
sample_nodes[[i]] <- remaining_node_numbers
}
}

write(unlist(sample_nodes), "thirty_nodes.txt")

for (j in 1:length(sample_nodes)){
for (i in 1:length(sample_nodes[[j]])){

full_extraction <- extract.clade(extant_tree, sample_nodes[[j]][[i]])
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

sample_thirty_nodes <- vector("list", length(species))
for (i in 1:length(species)){
sample_thirty_nodes[[i]] <- unique(species[[i]])
}

##########THREE_TAXON_TREES

two_node_index <- vector("list", length(sample_two_nodes))
for (j in 1:length(sample_two_nodes)){
for (i in 1:length(sample_two_nodes[[j]])){
two_node_index[[j]] <- append(two_node_index[[j]], which(extant_tree$tip.label == sample_two_nodes[[j]][[i]])) 
}
}

three_taxon_trees <- vector("list", length(two_node_index))
for (i in 1:length(three_taxon_trees)){
three_taxon_trees[[i]] <- drop.tip(extant_tree, extant_tree$tip.label[-two_node_index[[i]]])
}

##########THIRTY_TAXON_TREES

thirty_node_index <- vector("list", length(sample_thirty_nodes))
for (j in 1:length(sample_thirty_nodes)){
for (i in 1:length(sample_thirty_nodes[[j]])){
thirty_node_index[[j]] <- append(thirty_node_index[[j]], which(extant_tree$tip.label == sample_thirty_nodes[[j]][[i]])) 
}
}

thirty_taxon_trees <- vector("list", length(thirty_node_index))
for (i in 1:length(thirty_taxon_trees)){
thirty_taxon_trees[[i]] <- drop.tip(extant_tree, extant_tree$tip.label[-thirty_node_index[[i]]])
}
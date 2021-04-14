#' generates three taxon trees and three taxon sequence data
#' @export

GenerateThreeTaxonData <- function(root_age_tree, bl_tree, alignment_dir, exhaustive){  

dir.create("three_taxon_data")

####################################################################
####################################################################
###GET_KEY_INFO#####################################################
####################################################################
####################################################################

###GET_RELEVENT_SUBSET_TIPS#########################################

half_one <- extract.clade(root_age_tree, length(root_age_tree$tip.label)+2)$tip.label
half_two <- drop.tip(root_age_tree, half_one)$tip.label

print("got subset tips for root")

###TIP_NODE_NUMBERS#################################################

half_one_tip_node_numbers <- which(root_age_tree$tip.label %in% half_one)
half_two_tip_node_numbers <- which(root_age_tree$tip.label %in% half_two)

print("got tip node numbers for subset tips for root")

###EXHAUSTIVE_SAMPLING_COMBINATIONS_OF_INPUTS#######################

if (length(half_one_tip_node_numbers) > 1 & length(half_two_tip_node_numbers) > 1){ 
half_one_selections <- combn(half_one_tip_node_numbers, 2)
half_two_selections <- combn(half_two_tip_node_numbers, 2)
}

if (length(half_one_tip_node_numbers) > 1 & length(half_two_tip_node_numbers) == 1){ 
half_one_selections <- combn(half_one_tip_node_numbers, 2)
half_two_selections <- data.frame(half_two_tip_node_numbers)
}

if (length(half_one_tip_node_numbers) == 1 & length(half_two_tip_node_numbers) > 1){ 
half_one_selections <- data.frame(half_one_tip_node_numbers)
half_two_selections <- combn(half_two_tip_node_numbers, 2)
}

print("initiated exhaustive sampling combinations")

###################################################################
###################################################################
###GET_TIP_SAMPLES#################################################
###################################################################
###################################################################

###EXHAUSTIVE######################################################

if (exhaustive == "FULL"){

print("doing exhaustive subsetting")

if (nrow(half_one_selections) > 0){
half_one_tip_subsets <- vector("list", 0)
for (j in 1:ncol(half_one_selections)){
for (i in 1:length(half_two_tip_node_numbers)){
half_one_tip_subsets[[length(half_one_tip_subsets)+1]] <- c(root_age_tree$tip.label[half_one_selections[,j]], root_age_tree$tip.label[half_two_tip_node_numbers[[i]]])
}
}
}

if (nrow(half_two_selections) > 0){
half_two_tip_subsets <- vector("list", 0)
for (j in 1:ncol(half_two_selections)){
for (i in 1:length(half_one_tip_node_numbers)){
half_two_tip_subsets[[length(half_two_tip_subsets)+1]] <- c(root_age_tree$tip.label[half_two_selections[,j]], root_age_tree$tip.label[half_one_tip_node_numbers[[i]]])
}
}
}

}

###GENERATE_EXHAUSTIVE_INGROUP_SAMPLE#############################

if (exhaustive == "INGROUP"){

print("doing subsetting with ingroup exhaustive")

outgroup_from_half_one <- half_one_tip_node_numbers[[which(node.depth.edgelength(bl_tree)[half_one_tip_node_numbers] == sort(node.depth.edgelength(bl_tree)[half_one_tip_node_numbers])[[round((length(half_one_tip_node_numbers)/2), 0)+1]])]] 
outgroup_from_half_two <- half_two_tip_node_numbers[[which(node.depth.edgelength(bl_tree)[half_two_tip_node_numbers] == sort(node.depth.edgelength(bl_tree)[half_two_tip_node_numbers])[[round((length(half_two_tip_node_numbers)/2), 0)+1]])]] 

if (nrow(half_one_selections) > 0){
half_one_tip_subsets <- vector("list", 0)
for (i in 1:ncol(half_one_selections)){
half_one_tip_subsets[[length(half_one_tip_subsets)+1]] <- c(root_age_tree$tip.label[half_one_selections[,i]], root_age_tree$tip.label[outgroup_from_half_two])
}
}

if (nrow(half_two_selections) > 0){
half_two_tip_subsets <- vector("list", 0)
for (i in 1:ncol(half_two_selections)){
half_two_tip_subsets[[length(half_two_tip_subsets)+1]] <- c(root_age_tree$tip.label[half_two_selections[,i]], root_age_tree$tip.label[outgroup_from_half_one])
}
}

}

###GENERATE_EXHAUSTIVE_OUTGROUP_SAMPLE#############################

if (exhaustive == "OUTGROUP"){

print("doing subsetting with outgroup exhaustive")

half_one_tip_subsets <- vector("list", 0)
if (nrow(half_one_selections) > 1){
half_one_tree <- extract.clade(root_age_tree, findMRCA(root_age_tree, half_one, "node"))
for(j in 1:(length(half_one)-1)){
extraction <- extract.clade(half_one_tree, length(half_one) + j)
if (length(extraction$tip.label) > 2){
extraction_half_one <- sample(extract.clade(extraction, length(extraction$tip.label)+2)$tip.label, 1)
extraction_half_two <- sample(extraction$tip.label[-which(extraction$tip.label %in% extraction_half_one)], 1)
for (i in 1:length(half_two_tip_node_numbers)){
half_one_tip_subsets[[length(half_one_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[half_two_tip_node_numbers[[i]]]], extraction_half_one, extraction_half_two)
}
} else {
for (i in 1:length(half_two_tip_node_numbers)){
half_one_tip_subsets[[length(half_one_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[half_two_tip_node_numbers[[i]]]], extraction$tip.label[[1]], extraction$tip.label[[2]])
}
}
}
}


half_two_tip_subsets <- vector("list", 0)
if (nrow(half_two_selections) > 1){
half_two_tree <- extract.clade(root_age_tree, findMRCA(root_age_tree, half_two, "node"))
for(j in 1:(length(half_two)-1)){
extraction <- extract.clade(half_two_tree, length(half_two) + j)
if (length(extraction$tip.label) > 2){
extraction_half_one <- sample(extract.clade(extraction, length(extraction$tip.label)+2)$tip.label, 1)
extraction_half_two <- sample(extraction$tip.label[-which(extraction$tip.label %in% extraction_half_one)], 1)
for (i in 1:length(half_one_tip_node_numbers)){
half_two_tip_subsets[[length(half_two_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[half_one_tip_node_numbers[[i]]]], extraction_half_one, extraction_half_two)
}
} else {
for (i in 1:length(half_one_tip_node_numbers)){
half_two_tip_subsets[[length(half_two_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[half_one_tip_node_numbers[[i]]]], extraction$tip.label[[1]], extraction$tip.label[[2]])
}
}
}
}

}

###GENERATE_INEXHAUSTIVE_SAMPLING#################################

if (exhaustive == "NO"){

print("doing inexhaustive subsetting")

outgroup_from_half_one <- half_one_tip_node_numbers[[which(node.depth.edgelength(bl_tree)[half_one_tip_node_numbers] == sort(node.depth.edgelength(bl_tree)[half_one_tip_node_numbers])[[round((length(half_one_tip_node_numbers)/2), 0)+1]])]]
outgroup_from_half_two <- half_two_tip_node_numbers[[which(node.depth.edgelength(bl_tree)[half_two_tip_node_numbers] == sort(node.depth.edgelength(bl_tree)[half_two_tip_node_numbers])[[round((length(half_two_tip_node_numbers)/2), 0)+1]])]]

half_one_tip_subsets <- vector("list", 0)
if (nrow(half_one_selections) > 1){
half_one_tree <- extract.clade(root_age_tree, findMRCA(root_age_tree, half_one, "node"))
for(i in 1:(length(half_one_tree$tip.label)-1)){
extraction <- extract.clade(half_one_tree, length(half_one_tree$tip.label) + i)
if (length(extraction$tip.label) > 2){
extraction_half_one <- sample(extract.clade(extraction, length(extraction$tip.label)+2)$tip.label, 1)
extraction_half_two <- sample(extraction$tip.label[-which(extraction$tip.label %in% extraction_half_one)], 1)
half_one_tip_subsets[[length(half_one_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[outgroup_from_half_two]], extraction_half_one, extraction_half_two)
} else {
half_one_tip_subsets[[length(half_one_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[outgroup_from_half_two]], extraction$tip.label[[1]], extraction$tip.label[[2]])
}
}
}


half_two_tip_subsets <- vector("list", 0)
if (nrow(half_two_selections) > 1){
half_two_tree <- extract.clade(root_age_tree, findMRCA(root_age_tree, half_two, "node"))
for(i in 1:(length(half_two_tree$tip.label)-1)){
extraction <- extract.clade(half_two_tree, length(half_two_tree$tip.label) + i)
if (length(extraction$tip.label) > 2){
extraction_half_one <- sample(extract.clade(extraction, length(extraction$tip.label)+2)$tip.label, 1)
extraction_half_two <- sample(extraction$tip.label[-which(extraction$tip.label %in% extraction_half_one)], 1)
half_two_tip_subsets[[length(half_two_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[outgroup_from_half_one]], extraction_half_one, extraction_half_two)
} else {
half_two_tip_subsets[[length(half_two_tip_subsets) + 1]] <- c(root_age_tree$tip.label[[outgroup_from_half_one]], extraction$tip.label[[1]], extraction$tip.label[[2]])
}
}
}

}

##################################################################
##################################################################
###WRITE_SEQUENCE_SUBSETS_FOR_EACH_GENE_TO_FILE###################
##################################################################
##################################################################

`%notin%` <- Negate(`%in%`)
all_tip_subsets <- c(half_one_tip_subsets, half_two_tip_subsets)

for (j in 1:length(all_tip_subsets)){
print(paste("sorting three taxon tree", j, sep=" "))
setwd("three_taxon_data")
dir.create(paste("three_taxon_", j, sep=""))
setwd(paste("../", alignment_dir, sep=""))
for (i in 1:length(list.files(pattern = "\\.fasta$"))){
alignment <- read.phyDat(list.files(pattern = "\\.fasta$")[[i]], format = "fasta")
gene_name <- substr(list.files(pattern = "\\.fasta$")[[i]], 1, which(unlist(strsplit(list.files(pattern = "\\.fasta$")[[i]], "")) == ".")-1) 
setwd(paste("../", paste("three_taxon_data/three_taxon_", j, sep=""),sep=""))
if (length(which(all_tip_subsets[[j]] %in% names(alignment))) == 3){
write.phyDat(subset(alignment, all_tip_subsets[[j]]), paste(gene_name, ".nexus", sep=""), format = "nexus")
}
setwd(paste("../../",alignment_dir,sep=""))
}
setwd(paste("../", paste("three_taxon_data/three_taxon_", j, sep=""),sep=""))
write.tree(drop.tip(root_age_tree, which(root_age_tree$tip.label %notin% all_tip_subsets[[j]])), "three_taxon.tre")
setwd("../../")
}

print("extracted_sequence_subsets")

##################################################################
##################################################################
###CONCATENATE_SEQUENCES##########################################
##################################################################
##################################################################

for (j in 1:length(all_tip_subsets)){
setwd(paste("three_taxon_data/three_taxon_", j, sep=""))
phyDat_list <- vector("list", length(list.files(pattern = "\\.nexus$")))
for (i in 1:length(list.files(pattern = "\\.nexus$"))){
phyDat_list[[i]] <- read.phyDat(list.files(pattern = "\\.nexus$")[[i]], format = "nexus")
}
phyDat_concat <- phyDat_list[[1]]
for (i in 2:length(phyDat_list)){
phyDat_concat <- append(phyDat_concat, phyDat_list[[i]])
}
write.phyDat(phyDat_concat, "concatenated.nexus", format = "nexus")
setwd("../..")
}

print("complete")

##################################################################
##################################################################
##################################################################
##################################################################
##################################################################

}

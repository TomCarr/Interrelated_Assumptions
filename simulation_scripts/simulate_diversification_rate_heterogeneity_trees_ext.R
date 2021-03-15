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
lambda <- c(21, 21)
mu <- c(0, 18)
frac <- c(1, 1)
times <- c(0, 0.1)
n_reps <- 1
K <- 0

##########

tip_vector <- seq(1, n_tips, 1)
node_numbers <- seq(length(tip_vector)+1, (length(tip_vector)+2) + (length(tip_vector)-3), 1) 

##########

test_age <- 0
while (test_age <= 0.4){
entire_tree <- sim.rateshift.taxa(n_tips, n_reps, lambda, mu, frac, times, complete = TRUE, K=0)
entire_tree <- entire_tree[[1]]
test_age <- max(node.depth.edgelength(drop.tip(entire_tree, getExtinct(entire_tree))))
}

write.tree(entire_tree, "1_entire.tre")
write.tree(drop.tip(entire_tree, getExtinct(entire_tree)), "1.tre")

##########


# BasicTimeEstimates

This package enables node age estimation in phylogenies using a series of three taxon phylogenies. The three taxon phylogenies are sampled from the main phylogeny, such that every node in the main phylogeny is sampled.
The purpose of these analyses is to enable estimation of divergence times in a manner that is not influenced by the tree prior, which is often biologically unrealistic.
Divergence time estimates derived from this method will therefore primarily be sensitive to the molecular clock model.

Requires: phangorn, devtools, phytools, phylobase

There are currently 2 functions:

## Generate data subsets for inference of three taxon trees
GenerateThreeTaxonData(root_age_tree, bl_tree, alignment_dir, exhaustive)\
\
### Arguments:\
root_age_tree: an ultrametric species tree (with all taxa), with the root constrained to a single point calibration\
bl_tree: a species tree (with all taxa) with molecular branch lengths\
alignment_dir: directory of individual gene alignments (for all taxa)\
exhaustive: sampling strategy used to generate the three taxon data subsets. With "FULL", every single combination of tips is sampled. With "INGROUP", every single combination of tips in the ingroup is sampled, but the same outgroup is used. With "OUTGROUP" every single combination of outgroups is sampled, but the ingroup sampling is only designed such that every node is sampled. With "NO" sampling is designed such that every node in the overall species tree is sampled just once.\
### Outputs:\
three taxon alignments, and start trees\

## Compile age estimates from a large number of three taxon trees
CompileAgeEstimates(n_three_taxon_trees, species_tree, outputs_dir)\
\
### Arguments:\
number_three_taxon_trees: number of three taxon trees estimated\
species_tree: the species tree\
outputs_dir: output directory of the estimated three taxon tree\




# BasicTimeEstimates 

This r package enables node age estimation in phylogenies using a series of three taxon phylogenies. The three taxon phylogenies are derived from the main phylogeny, such that every node in the main phylogeny is sampled.
The purpose of these analyses is to enable estimation of divergence times in a manner that is not influenced by the tree prior, which is often biologically unrealistic.
Divergence time estimates derived from this method will therefore be primarily sensitive to the molecular clock model.

Fossil calibrations at internal nodes cannot currently be implemented with this method. Its primary purpose is the case where there are very few calibrations and age estimates are seisitive to the tree prior. However, this functionality will be added shortly.

Requires: phangorn, devtools, phytools, phylobase

There are currently 2 functions:

## Generate data subsets for inference of three taxon trees
GenerateThreeTaxonData(root_age_tree, bl_tree, alignment_dir, exhaustive)\

### Arguments:
root_age_tree: an ultrametric species tree (with all taxa), with the root constrained to a single point calibration. (Phylo object)\
bl_tree: a species tree (with all taxa) with molecular branch lengths (Phylo object)\
alignment_dir: directory of individual gene alignments (for all taxa). Alignments must be in fasta format\
exhaustive: sampling strategy used to generate the three taxon data subsets. With "FULL", every single combination of tips is sampled. With "INGROUP", every single combination of tips in the ingroup is sampled, but the same outgroup is used. With "OUTGROUP" every single combination of outgroups is sampled, but the ingroup sampling is only designed such that each node in the main phylogeny is sampled once. With "NO" sampling is designed such that every node in the overall species tree is sampled just once.
### Outputs:
three taxon alignments (in nexus format), and start trees for each three taxon alignment

## Compile age estimates from a large number of three taxon trees
CompileAgeEstimates(n_three_taxon_trees, species_tree, outputs_dir)

### Arguments:
number_three_taxon_trees: number of three taxon trees estimated\
species_tree: the species tree\
outputs_dir: output directory of the estimated three taxon trees
### Outpupts:
list of clade age estimates

# RevBayes scripts
Example scripts for estimating divergence times across many three taxon trees






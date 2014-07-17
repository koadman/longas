#!/bin/bash

FASTTREE="scripts/PhyloSift/bin/FastTree"
PDA="scripts/PhyloSift/bin/pda"

# build trees from the alignments
cat data/full_length/all_clean_try5.V4.full.aligned.fasta | $FASTTREE -nt -gtr > output/full_length/long_cmalign.tre
cat data/full_length/all_clean_try5.V4.aligned.fasta | $FASTTREE -nt -gtr > output/full_length/V4_cmalign.tre
cat data/full_length/all_clean_try5.V14.aligned.fasta | $FASTTREE -nt -gtr > output/full_length/V14_cmalign.tre

# calculate the phylogenetic diversity in the tree (just summing branch length of all taxa)
$PDA -k 999999999 output/full_length/long_cmalign.tre
$PDA -k 999999999 output/full_length/V4_cmalign.tre
$PDA -k 999999999 output/full_length/V14_cmalign.tre

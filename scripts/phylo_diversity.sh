#!/bin/bash

MUSCLE="~/software/muscle/muscle"
FASTTREE="~/git/PhyloSift/bin/FastTree"
PDA="~/git/PhyloSift/bin/pda"

# make the multiple alignments
$MUSCLE -in data/full_length/long_take_1_2_1000bp.fasta -out output/full_length/long_take_1_2_1000bp_aln.fa
$MUSCLE -in data/full_length/V4_take_1_2_1000bp.fasta -out output/full_length/V4_take_1_2_1000bp_aln.fa
bzip2 output/full_length/long_take_1_2_1000bp_aln.fa
bzip2 output/full_length/V4_take_1_2_1000bp_aln.fa

# build trees from the alignments
bzcat output/full_length/long_take_1_2_1000bp_aln.fa.bz2 | $FASTTREE -nt -gtr > output/full_length/long_take_1_2_1000bp.tre
bzcat output/full_length/V4_take_1_2_1000bp_aln.fa.bz2 | $FASTTREE -nt -gtr > output/full_length/V4_take_1_2_1000bp.tre

# calculate the phylogenetic diversity in the tree (just summing branch length of all taxa)
$PDA -k 999999999 output/full_length/long_take_1_2_1000bp.tre
$PDA -k 999999999 output/full_length/V4_take_1_2_1000bp.tre

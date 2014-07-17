#!/usr/bin/env python
#
# extract the V4 region from an infernal alignment
# assumes non-model columns have been masked and all sequences are on a single line 
#
import sys
from string import maketrans 

cmalign_infile = open(sys.argv[1], 'r')
v4_outfile = open(sys.argv[2], 'w')
v4_full_outfile = open(sys.argv[3], 'w')
v14_outfile = open(sys.argv[4], 'w')

line_i = -1
name_line = ""
for line in cmalign_infile:
	line_i+=1
	if line_i % 2 == 0:
		if line[0] != ">":
			sys.stderr.write("Error parsing " + sys.argv[1] + ".\nIt appears that alignments may span multiple lines\n")
		line = line.translate(maketrans(':','.'))
		name_line = line
	else:
		# TODO: alignment of V4 to infernal model suggests position 554 for start, 809 for end
		subaln = line[554:809]
		full_subaln = line[40:1410]
		v14_subaln = line[40:809]
		gapct = 0
		for i in range(len(subaln)):
			if(subaln[i] == '-'):
				gapct+=1
		if(gapct > len(subaln) / 2):
			continue	# skip seqs that are > 50% gap in the V4

		gapct = 0
		for i in range(len(v14_subaln)):
			if(v14_subaln[i] == '-'):
				gapct+=1
		if(gapct > len(v14_subaln) / 2):
			continue	# skip seqs that are > 50% gap in the V1-4


		gapct = 0
		for i in range(len(full_subaln)):
			if(full_subaln[i] == '-'):
				gapct+=1
		if(gapct > len(full_subaln) / 20):
			continue	# skip seqs that are > 5% gap overall

		v4_outfile.write(name_line + subaln + "\n")
		v14_outfile.write(name_line + v14_subaln + "\n")
		v4_full_outfile.write(name_line + full_subaln + "\n")


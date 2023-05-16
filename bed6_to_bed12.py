#!/usr/bin/env python3

# utility to convert sorted bed6 to bed12 - combines features based on shared name (bed column 4)

# 1=bed6 file input, 2=bed12 file output

import sys

bed12_dict={}

with open(sys.argv[1]) as bed6: # must be sorted bed6
	for line in bed6:
		line=line.rstrip().split('\t')
		if line[3] not in bed12_dict: 
			block_size=str(int(line[2])-int(line[1]))+','
			add=[line[0], line[1], line[2], line[3], line[4], line[5], line[1], line[2], '255,0,0', 1, block_size, '0,']
			bed12_dict[line[3]]=add
		else: # for genes with multiple exons
			change=bed12_dict[line[3]]
			change[2]=line[2] # change end coordinate
			change[7]=line[2] # change thick end coordinate
			change[9]+=1 # add second exon
			block_size=str(int(line[2])-int(line[1]))
			block_start=int(line[1])-int(change[1])
			change[10]=str(change[10])+str(block_size)+','
			change[11]=str(change[11])+str(block_start)+','
			bed12_dict[line[3]]=change

with open(sys.argv[2],'w') as outfile:
	for gene, bed_coords in bed12_dict.items():
		write_out=bed_coords
		write_out[9]=str(write_out[9])
		outfile.write('\t'.join(write_out)+'\n')

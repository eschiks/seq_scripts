#!/usr/bin/env python3

# inputs: 1=gtf file, 2=SJ.out.tab file from STAR (with unique and multi-mapped columns)
# prints out % of unannotated junctions 

import sys

# pulls transcript name from gtf description line

def find_tx_id(x):
    start = x.index('transcript_id') + len('transcript_id')+2
    start_to_end = x[start:]
    end = start_to_end.index(';')-1
    return start_to_end[:end]

# pulls introns from exon start and end coords for given transcript

def pull_intron_coords(exon_coords):
    result=[]
    if len(exon_coords)<2:
        return False
    else:
        for idx, x in enumerate(exon_coords):
            if idx==len(exon_coords)-1:
                continue
            else:
                start=x[1]+1
                end=exon_coords[idx+1][0]-1
                if start<=end:
                    result.append((start,end))
                else:
                    return False
        return result

# make dictionary with exon start and end coords for every transcript	

exon_dict={}
with open(sys.argv[1], 'r') as gtf: # gtf annotation
    for line in gtf:
        if line.startswith("#"):
            pass
        else:
            line=line.rstrip().split('\t')
            if line[2]=='exon':
                tx,chrom,exon_start,exon_end=find_tx_id(line[8]),line[0],int(line[3]),int(line[4])
                if (tx,chrom) not in exon_dict:
                    exon_dict[(tx,chrom)]=[(exon_start,exon_end)]
                else:
                    exon_dict[(tx,chrom)].append((exon_start,exon_end))

# parse exon dictionary to get intron coordinates - retain unique coordinates by adding to set

intron_set=set()
for gene in exon_dict:
    if not pull_intron_coords(sorted(exon_dict[gene])):
        continue
    else:
        tx=gene[0]
        chrom=gene[1]
        introns=(pull_intron_coords(sorted(exon_dict[gene])))
        for x in introns:
            intron_set.add((chrom,x[0],x[1]))

# go through junction file from STAR output

annotated=0
unannotated=0

with open(sys.argv[2], 'r') as junctions:
	for line in junctions:
		line=line.rstrip().split('\t')
		chrom,start,end,counts=line[0],int(line[1]),int(line[2]),int(line[6])
		if (chrom,start,end) in intron_set:
			annotated+=counts
		else:
			unannotated+=counts

total=unannotated+annotated
percent_unannotated=(float(unannotated)/float(total))*100
print(percent_unannotated)

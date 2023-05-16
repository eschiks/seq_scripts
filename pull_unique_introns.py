#!/usr/bin/env python3

# pulls intron coordinates from Wormbase gtf file that do not overlap with any annotated exon and outputs a SAF file that can be used with featureCounts

# inputs: 1=gtf file, 2=output bed file name, 3=output saf file name #

import sys,pybedtools,os

# pulls gene name from gtf description line

def find_wb_id(x):
    start = x.index('gene_id') + len('gene_id')+2
    start_to_end = x[start:]
    end = start_to_end.index(';')-1
    return start_to_end[:end]

# removes redundant exon coordinates

def remove_exon_overlaps(ranges):
    result = []
    current_start = -1
    current_stop = -1 
    for start, stop in sorted(ranges):
        if start > current_stop:
            result.append((start, stop))
            current_start, current_stop = start, stop
        else:
            result[-1] = (current_start, stop)
            current_stop = max(current_stop, stop)

    return result

# pulls intron coordinates from non-redundant exon coordinates

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

pos_exon_dict={}
neg_exon_dict={}

# pull all exon coordinates from gtf gene by gene

with open(sys.argv[1], 'r') as gtf:
    for line in gtf:
        if line.startswith("#"):
            pass
        else:
            line=line.rstrip().split('\t')
            if line[2]=='exon':
                wb,chrom,exon_start,exon_end,strand=find_wb_id(line[8]),line[0],int(line[3]),int(line[4]),line[6]
                if strand=='+':
                    if (wb,chrom,strand) not in pos_exon_dict:
                        pos_exon_dict[(wb,chrom,strand)]=[(exon_start,exon_end)]
                    else:
                        pos_exon_dict[(wb,chrom,strand)].append((exon_start,exon_end))
                elif strand=='-':
                    if (wb,chrom,strand) not in neg_exon_dict:
                        neg_exon_dict[(wb,chrom,strand)]=[(exon_start,exon_end)]
                    else:
                        neg_exon_dict[(wb,chrom,strand)].append((exon_start,exon_end))

# remove overlapping exon coordinates so each position is present only once

pos_non_overlap_exon_dict={}
neg_non_overlap_exon_dict={}

for gene in pos_exon_dict:
    pos_non_overlap_exon_dict[gene]=remove_exon_overlaps(pos_exon_dict[gene])

for gene in neg_exon_dict:
    neg_non_overlap_exon_dict[gene]=remove_exon_overlaps(neg_exon_dict[gene])
     
# pull intron coordinates from annotated exons (features that are always introns, never part of an exon)
        
pos_intron_dict={}
neg_intron_dict={}

for gene in pos_non_overlap_exon_dict:
    if not pull_intron_coords(pos_non_overlap_exon_dict[gene]):
        continue
    else:
        pos_intron_dict[gene]=pull_intron_coords(pos_non_overlap_exon_dict[gene])

for gene in neg_non_overlap_exon_dict:
    if not pull_intron_coords(neg_non_overlap_exon_dict[gene]):
        continue
    else:
        neg_intron_dict[gene]=pull_intron_coords(neg_non_overlap_exon_dict[gene])

# Write out temporary intron and exon bed files        

with open("tmp.bed",'w') as tmp:
    for gene in pos_intron_dict:
        name,chrom,strand=gene[0],gene[1],gene[2]
        for x in pos_intron_dict[gene]:
            start,end,score=str(x[0]),str(x[1]),str(60)
            tmp.write(chrom+'\t'+start+'\t'+end+'\t'+name+'\t'+score+'\t'+strand+'\n')
    for gene in neg_intron_dict:
        name,chrom,strand=gene[0],gene[1],gene[2]
        for x in neg_intron_dict[gene]:
            start,end,score=str(x[0]),str(x[1]),str(60)
            tmp.write(chrom+'\t'+start+'\t'+end+'\t'+name+'\t'+score+'\t'+strand+'\n')

with open("tmp2.bed",'w') as tmp2:
    for gene in pos_exon_dict:
        name,chrom,strand=gene[0],gene[1],gene[2]
        for x in pos_non_overlap_exon_dict[gene]:
            start,end,score=str(x[0]),str(x[1]),str(60)
            tmp2.write(chrom+'\t'+start+'\t'+end+'\t'+name+'\t'+score+'\t'+strand+'\n')
    for gene in neg_exon_dict:
        name,chrom,strand=gene[0],gene[1],gene[2]
        for x in neg_non_overlap_exon_dict[gene]:
            start,end,score=str(x[0]),str(x[1]),str(60)
            tmp2.write(chrom+'\t'+start+'\t'+end+'\t'+name+'\t'+score+'\t'+strand+'\n')

# Create a bed file that removes intronic sequences that overlap with an exonic sequence. Write it out as a bed file and a .saf file (for feature counts)

a = pybedtools.BedTool('tmp.bed') # introns file
b = pybedtools.BedTool('tmp2.bed') # exons file
intersect=a.intersect(b, v=True,s=True,output=sys.argv[2]) # outputs intron bed file

os.remove('tmp.bed')
os.remove('tmp2.bed')

# convert bed file into saf format 

with open(sys.argv[3],'w') as saf:
    saf.write('GeneID'+'\t'+'Chr'+'\t'+'Start'+'\t'+'End'+'\t'+'Strand'+'\n')
    with open(sys.argv[2]) as intersect:
        for line in intersect:
            line=line.rstrip().split('\t')
            chrom,start,end,name,strand=line[0],line[1],line[2],line[3],line[5]
            out=[name,chrom,str(start),str(end),strand]
            saf.write('\t'.join(out)+'\n')

#!/usr/bin/env python3

# pulls intergenic regions from Wormbase gtf file and outputs SAF file with coordinates to be used for featureCounts

# inputs: 1=chrom sizes file, 2=gtf file, 3=output saf file name #

import sys

def find_wb_id(x):
    start=x.index('gene_id')+len('gene_id')+2
    start_to_end=x[start:]
    end=start_to_end.index(';')-1
    return start_to_end[:end]

def assign_chrom_sizes(chrom_sizes_file):
    with open(chrom_sizes_file,'r') as chrom_sizes:
        chrom_size_dict={}
        for line in chrom_sizes:
            line=line.rstrip().split('\t')
            chrom_size_dict[line[0]]=line[1]
        return chrom_size_dict
    
def pull_left_and_right_exon(exon_list):
    sorted_list=sorted(exon_list)
    left=sorted_list[0]
    right=sorted_list[-1]
    chrom=left[0]
    left_coord=left[1]
    right_coord=right[2]
    return chrom,left_coord,right_coord

def remove_overlaps(ranges):
    result=[]
    current_start=-1
    current_end=-1
    for start,stop in sorted(ranges):
        if stop>=current_end:
            result.append((start,stop))
            current_start=start
            current_end=stop
    return result

# read in chromosome lengths and assign to variables

chrom_sizes=assign_chrom_sizes(sys.argv[1])

gene_dict={}

# pull exons from every gene into gene dictionary

with open(sys.argv[2], 'r') as gtf:
    for line in gtf:
        if line.startswith("#"):
            pass
        else:
            line=line.rstrip().split('\t')
            if line[2]=='exon' or line[2]=='five_prime_utr' or line[2]=='three_prime_utr': # some genes the UTRs aren't included as part of exons?
                wb,chrom,exon_start,exon_end=find_wb_id(line[8]),line[0],int(line[3]),int(line[4])
                if wb not in gene_dict:
                    gene_dict[wb]=[(chrom,exon_start,exon_end)]
                else:
                    gene_dict[wb].append((chrom,exon_start,exon_end))

chrom_dict={}

# pull the leftmost and rightmost coordinate from every gene

for gene in gene_dict:
    chrom,left,right=pull_left_and_right_exon(gene_dict[gene])
    if chrom not in chrom_dict:
        chrom_dict[chrom]=[(left,right)]
    else:
        chrom_dict[chrom].append((left,right))
        
# remove overlapping coordinates (i.e. if a gene lies entirely within another gene)

non_redundant={}

for x in chrom_dict:
    non_redundant[x]=remove_overlaps(chrom_dict[x])       

# pull intergenic regions
    
ig_dict={}
        
for y in non_redundant:
    result=[]
    for idx, x in enumerate(sorted(non_redundant[y])):
        if idx==len(sorted(non_redundant[y]))-1: # if last exon on chromosome
            start=x[1]+1
            end=int(chrom_sizes[y])
            if start<end:
                result.append((start,end))
        else:
            next_gene=sorted(non_redundant[y])[idx+1]
            if idx==0: # if first exon on chromosome, start is always 1, end is start of exon 1
                start0=1
                if x[0]>2:
                    end0=x[0]-1
                    result.append((start0,end0))
            if x[1]<next_gene[0]: # if there is a gap between genes
                start=x[1]+1
                end=next_gene[0]-1
                if end-start>1:
                    result.append((start,end))
    ig_dict[y]=result


with open(sys.argv[3],'w') as saf: 
    num=1 # give arbitrary names for featureCounts
    saf.write('GeneID'+'\t'+'Chr'+'\t'+'Start'+'\t'+'End'+'\t'+'Strand'+'\n')
    for x in ig_dict:
        chrom=x
        for y in ig_dict[chrom]:
            start,end=str(y[0]),str(y[1])
            strand='+' # this is arbitrary
            name='intergenic_region_'+str(num)
            num+=1
            saf.write(name+'\t'+chrom+'\t'+start+'\t'+end+'\t'+strand+'\n')

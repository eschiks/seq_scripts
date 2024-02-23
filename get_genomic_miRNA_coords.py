#!/usr/bin/env python3

# this script pulls all 3'UTR sequences from a GTF file and searches them for miRNA 7mer seeds from TargetScan 

# 1=sorted gtf file input, 2=genome fasta file input, 3=miRNA seed seq file (obtained here: https://www.targetscan.org/cgi-bin/targetscan/mirna_families.cgi?db=worm_52, keep first 2 columns) # 4=ouput bed file name

import sys
import pybedtools
import os

# make bed with 3'UTR coordinates from gtf annotated 3'UTRs

def find_tx_id(x):
	start = x.index('transcript_id') + len('transcript_id')+2
	start_to_end = x[start:]
	end = start_to_end.index(';')-1
	return start_to_end[:end]

def rev_complement(seq):
	rev_dict={'A':'t','U':'a','G':'c','C':'g'}
	for key in rev_dict:
		seq=seq.replace(key,rev_dict[key])
	rev_com=seq[::-1]
	return rev_com.upper()

# pull 3'UTR coordinates from gtf --> write tmp bed file

with open(sys.argv[1]) as gtf:
	with open('tmp.bed','w') as utr_bed:
		for line in gtf:
			if line.startswith('#'):
				continue
			else:
				line=line.strip().split('\t')
				if line[2]=='three_prime_utr':
					start=str(int(line[3])-1)
					if int(line[4])>int(line[3]):
						add=[line[0],start,line[4],find_tx_id(line[8]),'60',line[6]]
						utr_bed.write('\t'.join(add)+'\n')
# convert bed to tmp bed12 

bed12_dict={}

with open('tmp.bed') as bed6: # must be sorted bed6
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

utr_coords={}
			
with open('tmp.bed12','w') as outfile:
	for gene, bed_coords in bed12_dict.items():
		write_out=bed_coords
		write_out[9]=str(write_out[9])
		utr_coords[gene]=write_out # added
		outfile.write('\t'.join(write_out)+'\n')

a=pybedtools.BedTool('tmp.bed12')
fasta=pybedtools.example_filename(sys.argv[2])
a.sequence(fi=fasta,s=True,split=True,nameOnly=True).save_seqs('tmp.fa')

# make dictionary of reverse complement of seed sequences to search for in 3'UTRs

seven_mers={}
seed_list=[]

with open(sys.argv[3]) as mirs:
	for line in mirs:
		line=line.strip().split()
		mir,rev_string=line[0],rev_complement(line[1])
		seven_mers[rev_string]=mir
		seed_list.append(rev_string)

# create a dictionary of 3'UTR sequences separated by strand

utr_seqs_minus={}
utr_seqs_plus={}
with open('tmp.fa') as fasta:
	for line in fasta:
		line=line.strip()
		if line.startswith('>'):
			name=line[1:line.find('(')]
			strand=line[line.find('(')+1:line.find(')')]
		else:
			seq=line
			if strand=='-':
				utr_seqs_minus[name]=seq
			else:
				utr_seqs_plus[name]=seq

# get index of location of miRNA sites within 3'UTR sequences              

found_dict={}
for utr in utr_seqs_plus:
	seq=utr_seqs_plus[utr]
	for seed in seed_list:
		res = [i for i in range(len(seq)) if seq.startswith(seed, i)]
		if res:
			if seed not in found_dict:
				found_dict[seed]=[[utr,res]]
			else:
				found_dict[seed].append([utr,res])

for utr in utr_seqs_minus:
	seq=utr_seqs_minus[utr]
	for seed in seed_list:
		res = [i for i in range(len(seq)) if seq.startswith(seed, i)]
		if res:
			res_minus=[len(seq)-i-7 for i in res] # changed from -7
			if seed not in found_dict:
				found_dict[seed]=[[utr,res_minus]]
			else:
				found_dict[seed].append([utr,res_minus])

# get info from 3'UTRs to include in final bed output                

#utr_coords={}
#with open('tmp.bed12') as bed:
	#for line in bed:
		#line=line.strip().split()
		#utr_coords[line[3]]=line

# convert locations of 7mers relative to 3'UTR starts to genomic coordinates, accounting for the fact that some seeds may span exon-exon junctions
        
bed=[]
for seed in found_dict:
	mir=seven_mers[seed]
	sites=found_dict[seed]
	for utr in sites:
		name=utr[0]
		seed_coords=utr[1]
		for site in seed_coords:
			utr_start=utr_coords[name][6]
			if utr_coords[name][-1]=='0,':
				seed_start=str(int(utr_start)+int(site))
				seed_end=str(int(seed_start)+7)
				if utr_coords[name][0]=='MtDNA':
					chrom,strand='M',utr_coords[name][5]
				else:
					chrom,strand=utr_coords[name][0],utr_coords[name][5]
				add=[chrom,seed_start,seed_end,seed,'60',strand]
				if add not in bed:
					bed.append(add) # works!
			else:
				sizes=utr_coords[name][10].split(',')
				starts=utr_coords[name][11].split(',')
				for idx,size in enumerate(sizes):
					size=int(size)
					if site>size:
						site=site-size
					else:
						if utr_coords[name][0]=='MtDNA':
							chrom,strand='M',utr_coords[name][5]
						else:
							chrom,strand=utr_coords[name][0],utr_coords[name][5]
						seed_start=int(utr_coords[name][1])+int(starts[idx])+site
						end_site=site+7
						if end_site <= size: # if seed site ends in the same exon as it started
							seed_end=seed_start+7
							add=[chrom,str(seed_start),str(seed_end),seed,'60',strand]
							if add not in bed:
								bed.append(add)
						else: # if seed starts in one exon and ends in another
							seed_end=int(utr_coords[name][1])+int(starts[idx])+int(size)
							remaining=7-(seed_end-seed_start)
							seed_start_2=int(utr_coords[name][1])+int(starts[idx+1])
							seed_end_2=seed_start_2+remaining
							add,add2=[chrom,str(seed_start),str(seed_end),seed,str(remaining),strand],[chrom,str(seed_start_2),str(seed_end_2),seed,'60',strand]
							if add not in bed and add2 not in bed:
								bed.append(add)
								bed.append(add2)
						break 

# write out the genomic coordinates of 7mer seeds to a bed6 file

with open(sys.argv[4], 'w') as outfile:
	for x in bed:
		x[3]=seven_mers[x[3]]
		outfile.write('\t'.join(x)+'\n')
		
os.remove('tmp.bed')
os.remove('tmp.bed12')
os.remove('tmp.fa')

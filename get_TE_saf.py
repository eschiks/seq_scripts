#!/usr/bin/env python3

# pulls transposable elements from UCSC rmsk file and outputs SAF file for featureCounts

filter_list=['DNA','DNA?','LINE','LTR','RC','SINE'] # filter list from Garrigues et al. 2022 - all types of transposable elements

with open('rmsk_TEs.saf','w') as saf:
	saf.write('GeneID'+'\t'+'Chr'+'\t'+'Start'+'\t'+'End'+'\t'+'Strand'+'\n')
	with open('rmsk.txt', 'r') as rmsk:
		for line in rmsk:
			line=line.rstrip().split('\t')
			chrom,start,end,name,strand,repclass=line[5].replace('chr',''),line[6],line[7],line[10],line[9],line[11]
			if any(TE in repclass for TE in filter_list):
				out=[name,chrom,str(start),str(end),strand]
				saf.write('\t'.join(out)+'\n')
		
		

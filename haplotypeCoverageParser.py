#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
Script to parse the coverage of linked SNPs within genes (gene-wise haplotypes).

Input: 
Linkage-filtered gene-wise haplotypes for each gene (output from haplotype_parser_v3.py script).
File containing paths to .bam files (one file per line).

Output:
A table that can be imported to R, for example, containing the summarized counts (read coverage in .bam files) of "A" and "B" haplotypes per gene for each sample. 
'''

def get_ref(gene, pos):
	from Bio import SeqIO
    #ref = open('/home/jpverta/software/bowtie-1.0.0/indexes/GCAT_WS-3.3.cluseq.fa','r')
	#ref = open('/Users/peckaj/Documents/Koulujutut/RNAseq/GCAT_WS-3.3.cluseq.fa','r')
    ref = open('/Volumes/LaCie/Data/Katak/Megagametophyte/expression/GCAT_WS-3.3.cluseq.fa','r')
    for record in SeqIO.parse(ref,'fasta'):
		if record.id == gene:
			return record.seq[int(pos)-1] 
		else:
			pass


def get_genes(infile):
	with open(infile,'r') as file:
		for line in file:
			yield line.split('\t')[0], []



def haplo_dict(haplofile,max_SNP=True,max_distance=False):
	import re,ast
	haplo_dict = dict(get_genes(haplofile))
	if max_SNP == True and max_distance == False:
		for haploline in open(haplofile,'r').read().split('\n'):
			gene = haploline.split('\t')[0]
			if len(gene) == 0:
				break
			elif len(haplo_dict[gene]) == 0: # if a dictioanry entry doesnt exist for the gene
				haplo_dict[gene] = (
									[re.sub('([\[\] ])','',x) for x in haploline.split('\t')[1].split(',')],
									[re.sub('([\[\] ])','',x) for x in haploline.split('\t')[2].split(',')],
									[re.sub('([\'\[\] ])','',x) for x in haploline.split('\t')[3].split(',')],
									[ast.literal_eval('{'+x+'}') for x in re.findall(r'\{(.+?)\}',haploline.split('\t')[4])],
									[ast.literal_eval('{'+x+'}') for x in re.findall(r'\{(.+?)\}',haploline.split('\t')[5])],
									ast.literal_eval(haploline.split('\t')[6]),
									ast.literal_eval(haploline.split('\t')[7]))
			elif len(haplo_dict[gene]) > 0 and len(str(haploline.split('\t')[2]).split(',')) > len(haplo_dict[gene][2].split(',')): # if an entry exists but the new line has more SNPs in the haplotype
				haplo_dict[gene] = (
									[re.sub('([\[\] ])','',x) for x in haploline.split('\t')[1].split(',')],
									[re.sub('([\[\] ])','',x) for x in haploline.split('\t')[2].split(',')],
									[re.sub('([\'\[\] ])','',x) for x in haploline.split('\t')[3].split(',')],
									[ast.literal_eval('{'+x+'}') for x in re.findall(r'\{(.+?)\}',haploline.split('\t')[4])],
									[ast.literal_eval('{'+x+'}') for x in re.findall(r'\{(.+?)\}',haploline.split('\t')[5])],
									ast.literal_eval(haploline.split('\t')[6]),
									ast.literal_eval(haploline.split('\t')[7]))
			else: # if the entry exists and the line has less SNPs than the existing haplotype
				pass
	else:
		print 'in haplo_dict: you have to chose between max number of SNPs or max distance between SNPs in a haplotype'
	return haplo_dict

def haploFilter(haplodict,haplofile): # filter for SNPs that are separated with at least 100 bp, otherwise counts are not independent because single read can span over two SNPs
	filtered_haplodict = dict(get_genes(haplofile))
	for gene in haplodict.keys():
		filtered_haplodict[gene] = [[],[],[],[],[],{},{}]
		posInd=[]
		posSet=[]
		for SNP in haplodict[gene][1]:
			ind = haplodict[gene][1].index(SNP)
			if all([float(SNP)-float(y) > 100 for y in posSet]): # if the SNP positions already in posSet are separated with 100 bp
				posInd.append(ind) 
				posSet.append(SNP)
		filtered_haplodict[gene][0] = [haplodict[gene][0][i] for i in posInd]
		filtered_haplodict[gene][1] = [haplodict[gene][1][i] for i in posInd]
		filtered_haplodict[gene][2] = [haplodict[gene][2][i] for i in posInd]
		filtered_haplodict[gene][3] = [haplodict[gene][3][i] for i in posInd]
		filtered_haplodict[gene][4] = [haplodict[gene][4][i] for i in posInd]
		for s in posSet:
			filtered_haplodict[gene][5].update({int(s):haplodict[gene][5][int(s)]})
			filtered_haplodict[gene][6].update({int(s):haplodict[gene][6][int(s)]})
	return filtered_haplodict


def haplotypeCoverage(haplofile,bamfiles,outfile,alt_os_file='alt_observations_closeDrop.haplotypes',ref_os_file='ref_observations_closeDrop.haplotypes'):
	import re
	import numpy
	raw_haplodict = haplo_dict(haplofile)
	haplodict = haploFilter(raw_haplodict,haplofile)
	open(outfile,'a').write('%s\t%s\t%s\t%s\t' % ('Gene','Positions','Alt','Ref'))
	bamlines = open(bamfiles,'r').read().split('\n')
	samples = [re.sub("[\D]","",x) for x in bamlines if len(x) > 0]
	alt_os = dict.fromkeys(samples)
	ref_os = dict.fromkeys(samples)
	for s in samples:
		open(outfile,'a').write(str(s)+'_altcov'+'\t')
	for s in samples:
		open(outfile,'a').write(str(s)+'_refcov'+'\t')
	open(outfile,'a').write('\n')
	for gene in haplodict.keys():
		alt_cov = dict.fromkeys(samples)
		ref_cov = dict.fromkeys(samples)
		for SNP in haplodict[gene][1]:
			ind =  haplodict[gene][1].index(SNP)
			for sample in samples:
				if sample in haplodict[gene][5][int(SNP)]:
					try:
						if alt_os[sample]:
							alt_os[sample] = alt_os[sample] + haplodict[gene][3][ind][sample]
						else:
							alt_os[sample] = haplodict[gene][3][ind][sample]
						if alt_cov[sample]:
							alt_cov[sample] = alt_cov[sample] + haplodict[gene][3][ind][sample]
						else:
							alt_cov[sample] = haplodict[gene][3][ind][sample]
					except: # if sample is in list of alternative samples but its coverage is not in the coverage list 
						if alt_os[sample]:
							alt_os[sample] = alt_os[sample] + 0
						else:
							alt_os[sample] = 0
						if alt_cov[sample]:
							alt_cov[sample] = alt_cov[sample] + 0
						else:
							alt_cov[sample] = 0
				elif sample in haplodict[gene][6][int(SNP)]:
					try:
						if ref_os[sample]:
							ref_os[sample] = ref_os[sample] + haplodict[gene][4][ind][sample]
						else:
							ref_os[sample] = haplodict[gene][4][ind][sample]
						if ref_cov[sample]:
							ref_cov[sample] = ref_cov[sample] + haplodict[gene][4][ind][sample]
						else:
							ref_cov[sample] = haplodict[gene][4][ind][sample]
					except: # if sample is in list of alternative samples but its coverage is not in the coverage list 
						if ref_os[sample]:
							ref_os[sample] = ref_os[sample] + 0
						else:
							ref_os[sample] = 0
						if ref_cov[sample]:
							ref_cov[sample] = ref_cov[sample] + 0
						else:
							ref_cov[sample] = 0
				else:
					pass
		open(outfile,'a').write('%s\t%s\t%s\t%s\t' % (gene, haplodict[gene][1],haplodict[gene][2],[get_ref(gene,x) for x in haplodict[gene][1]]))
		for s in samples:
			open(outfile, 'a').write(str(alt_cov[s])+'\t')
		for s in samples:
			open(outfile, 'a').write(str(ref_cov[s])+'\t')
		open(outfile, 'a').write('\n')
	for sample in samples:
		open(alt_os_file,'a').write(sample+'\t'+str(alt_os[sample])+'\n')
		open(ref_os_file,'a').write(sample+'\t'+str(ref_os[sample])+'\n')


import sys 

try:
	haplofile = sys.argv[1]
	bamfiles = sys.argv[2]
	outfile = sys.argv[3]
	alt_os_file = sys.argv[4]
	ref_os_file = sys.argv[5]
except:
	print 'Check command line arguments!'

haplotypeCoverage(haplofile,bamfiles,outfile,alt_os_file,ref_os_file)


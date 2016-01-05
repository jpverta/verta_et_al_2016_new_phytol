#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014

Input: 
Gene-wise genotypes for each gene (output from genotype_parser_v2.py script).

Output:
A table that can be imported to R, for example, containing the summarized counts (read coverage in .bam files) per SNPs per gene for each sample. 
'''



def get_ref(gene, pos):
	from Bio import SeqIO
	ref = open('/home/jpverta/software/bowtie-1.0.0/indexes/GCAT_WS-3.3.cluseq.fa','r')
	for record in SeqIO.parse(ref,'fasta'):
		if record.id == gene:
			return record.seq[int(pos)-1] 
		else:
			pass


def get_genes(infile):
	with open(infile,'r') as file:
		for line in file:
			yield line.split('\t')[0], ()



def snp_dict_parser(snpfile):
        import re,ast
        snp_dict = dict(get_genes(snpfile))
        for snpline in open(snpfile,'r').read().split('\n'):
        	gene = snpline.split('\t')[0]
                if len(gene) == 0:
                	break
                elif len(snp_dict[gene]) == 0: # if a dictioanry entry doesnt exist for the gene
					snp_dict[gene] = (
					(re.sub('([\[\] ])','',snpline.split('\t')[1])),
					(re.sub('([\[\] ])','',snpline.split('\t')[2])),
					(re.sub('([\'\[\] ])','',snpline.split('\t')[3])),
					(ast.literal_eval(snpline.split('\t')[4])),
					[re.sub('([A-Za-z._\'\" \[\]])','',str(y)) for y in snpline.split('\t')[5].split(',')])
                else:
                	pass
        return snp_dict




def snpCoverage(snpfile,outfile,alt_os_file,ref_os_file):
	import re
	import numpy
	snp_dict = snp_dict_parser(snpfile)
	open(outfile,'a').write('%s\t%s\t%s\t%s\t' % ('Gene','Positions','Alt','Ref'))
	bamlines = open('filtered_bam_files.txt','r').read().split('\n')
	samples = [re.sub("[\D]","",x) for x in bamlines if len(x) > 0]
	alt_os = dict.fromkeys(samples)
	ref_os = dict.fromkeys(samples)
	for s in samples:
		open(outfile,'a').write(str(s)+'_altcov'+'\t')
	for s in samples:
		open(outfile,'a').write(str(s)+'_refcov'+'\t')
	open(outfile,'a').write('\n')
	for gene in snp_dict.keys():
		alt_cov = dict.fromkeys(samples)
		ref_cov = dict.fromkeys(samples)
		for sample in samples:
			if sample in [re.sub('[A-Za-z.]','',x) for x in snp_dict[gene][4]]:
				if alt_os[sample]:
					alt_os[sample] = alt_os[sample] + snp_dict[gene][3][sample]
				else:
					alt_os[sample] = snp_dict[gene][3][sample]
				if alt_cov[sample]:
					alt_cov[sample] = alt_cov[sample] + snp_dict[gene][3][sample]
				else:
					alt_cov[sample] = snp_dict[gene][3][sample]
			else:
				if ref_os[sample]:
					ref_os[sample] = ref_os[sample] + snp_dict[gene][3][sample]
				else:
					ref_os[sample] = snp_dict[gene][3][sample]
				if ref_cov[sample]:
					ref_cov[sample] = ref_cov[sample] + snp_dict[gene][3][sample]
				else:
					ref_cov[sample] = snp_dict[gene][3][sample]
		open(outfile,'a').write('%s\t%s\t%s\t%s\t' % (gene, snp_dict[gene][0],snp_dict[gene][2],get_ref(gene,snp_dict[gene][0][0])))
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
	snpfile = sys.argv[1]
	outfile = sys.argv[2]
	alt_os_file = sys.argv[3]
	ref_os_file = sys.argv[4]
except:
	print 'Check command line arguments!'

snpCoverage(snpfile,outfile,alt_os_file,ref_os_file)



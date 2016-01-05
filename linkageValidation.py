#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
'''

def position_coverage(sample, chromosome, position, bam_file_list=bamfiles,vcf_type=0):
        import pysam
        coverage = 0
        file_list = open(bam_file_list,'r').read().split('\n')
        for bam in file_list: # for each file in list
                if sample in bam and len(bam) > 0:
                        try:
                                samfile = pysam.Samfile(bam,'rb') # read BAM file with pysam
                                sampileup = samfile.pileup(chromosome, int(position)-1,int(position))
                                for pileupcolumn in sampileup: # look above concerning cordinates
                                         if vcf_type == 0: # if SNP positions in vcf file correspond exactly to positions in BAM file
                                                if pileupcolumn.pos == int(position): # look above concerning cordinates
                                                        coverage = int(pileupcolumn.n)
                                         elif vcf_type == -1: # if SNP positions iin VCF file are shifted -1 relative to BAM file (SAMTOOLS mpileup)
                                                if pileupcolumn.pos == position-1: # look above concerning cordinates
                                                        coverage = int(pileupcolumn.n)
                                         else:
                                                print 'Incorrect vcf_type -option, should be 0 or -1'
                                                break
                        except:
                                print 'Warning in position_coverage: Could not open BAM file %s' % bam
                                pass
                else:
                        pass
        return coverage,sample


def listComparison(ls,sample_lists):
	for o_l in sample_lists:
		return [x in o_l for x in ls]


def linkage(snpfile,outfile,bamfiles):
        import re
        import yaml 
        snp = yaml.load(open(snpfile,'r').read())
        bamlines = open(bamfiles,'r').read().split('\n')
        samples = [re.sub("[\D]","",x) for x in bamlines if len(x) > 0]
        for gene in snp.keys():
			if len(snp[gene]) > 0:
				cov = []
				alt_samples = []
				pos = []
				for SNPindex in range(len(snp[gene][0])):
					if len(snp[gene][4][SNPindex]) >= 25 and len(snp[gene][4][SNPindex]) <= 41: # criteria for SNP frequency, now set to Mendelian
						pos.append(snp[gene][0][SNPindex])
						alt_samples.append([re.sub("[\D]","",x) for x in snp[gene][4][SNPindex]])
					else:
						pass
				obs = [listComparison(x,alt_samples) for x in alt_samples] # a list o Booleans expressing wehter a SNP is observed in all samples (under linkage) 
				linkage =  [[sum(x),len(x)] for x in obs] 
				''' a list for each SNP, giving the number of samples which share the SNP and the number of compared samples. if they are equal it means that 
				the SNP is observed in each sample. if the first field is zero it means that it is an alternative SNP under complete linkage (the other allele of 
				77111). deviations between the first and second fields indicate that there is no linkage over samples '''
				all_linked = []
				for x in linkage:
					if x[1]-x[0] <= x[1]-0.9*x[1] or x[1]-x[0] >= 0.9*x[1]:
						all_linked.append(1)					
					else:
						all_linked.append(0)
				open(outfile,'a').write(gene+'\t'+str(pos)+'\t'+str(all_linked)+'\n')

import sys 

try:
	snpfile = sys.argv[1] # output from vcf_parser.py (.yaml file containing raw SNP variants)
	bamfiles = sys.argv[2] # names of filtered .bam files, one sample per line in text document 
except:
	print 'Check command line arguments!'
 
linkage(snpfile,'linkageValidation.txt',bamfiles)

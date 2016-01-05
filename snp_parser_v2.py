#!/home/jpverta/bin/python2.7

''' 
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014

Script to parse SNPs from a file containing all observed SNPs (output from vcf_parser.py -script).

Input:
Yaml -file containing the SNP dictionary (vcf_parser.py output)
Number of samples.
Confidence interval for segregation frequencies IC +- N/2.
Minimum coverage of SNP position required per sample.
Minimum number of samples where the coverage is attained.
File containing paths to .vcf files, one path per line.
File containing paths to .bam files, one path per line.

Output is a table of SNPs that meet the quality and segregation criteria set above. Linkage between SNPs in single genes will not be considered. Sets of linked SNPs are identified with haplotype_parser_v3.py script. 			 		
'''

def position_type(chromosome, position, vcf_file_list):
	import vcf
	type = []
	file_list = open(vcf_file_list,'r').read().split('\n')
	for file in file_list: # for each file in list
		try:	
			reader = vcf.Reader(filename=file)
			for record in reader.fetch(chromosome,position,position+1): # look above concerning cordinates
				type.append(record.INFO['TYPE'][0])
		except:
			print 'Warning in position_type: Could not open VCF file %s' % file
			pass
	return type


def position_coverage(chromosome, position, samples, bam_file_list):
        import pysam
		import re
        coverage = dict.fromkeys(samples)
        file_list = open(bam_file_list,'r').read().split('\n')
        for sample in samples:
                for bam in file_list: # for each file in list
                        if re.sub('[\D]','',bam) == sample:
                                samfile = pysam.Samfile(bam,'rb') # read BAM file with pysam
                                sampileup = samfile.pileup(chromosome, int(position)-1,int(position))
                                for pileupcolumn in sampileup: 
                                        if pileupcolumn.pos == position: 
                                                coverage[sample] = pileupcolumn.n
                if not coverage[sample]:
                        coverage[sample] = 0
        return coverage



def quality_filter(gene, pos, vcf_files):
        import vcf
        quality = []
        for vcf_file in open(vcf_files,'r').read().split('\n'):
                if len(vcf_file) > 1:
                        vcf_reader = vcf.Reader(filename=vcf_file)
                        try:
                                vcf_record = vcf_reader.fetch(gene, pos)
                                quality.append(vcf_record.INFO['MQM'][0])
                        except:
                                pass
        print 'quality ' + str(quality)
        return quality


def get_linkage(gene,linkageFile='linkageValidation.txt'):
        import re
        for line in open(linkageFile,'r').read().split('\n'):
                if line.split('\t')[0] == gene:
			linkage = dict(zip([re.sub("[\D]",'',x) for x in line.split('\t')[1].split(',')], [re.sub("[\D]",'',x) for x in line.split('\t')[2].split(',')]))
                	return linkage         




def snp_parser(in_yaml,outfile_parsed,vcf_files,bam_file_list,N,IC,min_coverage,min_samples,filter_complex=True):
	print in_yaml,outfile_parsed,vcf_files,bam_file_list,N,IC,min_coverage,min_samples
	import yaml
	import re
	open(outfile_parsed,'a')
	bases = []
	vcf = yaml.load(open(in_yaml,'r').read())
	samples =  [re.sub('[\D]','',bam) for bam in open(bam_file_list,'r').read().split('\n') if len(re.sub('[\D]','',bam)) > 1]
	for key in vcf.keys(): # for all chromosomes (transcripts)
		linkage = get_linkage(key)
		if len(vcf[key]) > 1 and len(linkage.keys()) >= 1 and len(linkage.keys()[0]) > 1 and (sum([x == '1' for x in linkage.values()]) >= 0.9*len(linkage.values())): # if it has SNPs and 90% the SNPs are under linkage
			for SNP in [vcf[key][0].index(int(SNP)) for SNP in linkage.keys() if linkage[SNP] == '1']:
				if vcf[key][3][SNP][0] in range(int(N)/2-int(IC) ,int(N)/2+int(IC)) and sum([qual >= 30 for qual in quality_filter(key,vcf[key][0][SNP],vcf_files)]) >= 0.9*len(quality_filter(key,vcf[key][0][SNP],vcf_files)):
					if len(vcf[key][3][SNP]) == 1: # if SNP not multiallelic
						if filter_complex == True: # if you have not defined non-parsed outfile, then check that you have +coverage (default 5) in each sampe at alt OR ref allele and all variants are SNPs '''
							coverage = position_coverage(key, vcf[key][0][SNP], samples, bam_file_list) # record the coverage of the position in all samples (in the directory)
							if sum([cov >= int(min_coverage) for cov in coverage.values()]) >= int(min_samples): # if the position is covered by min_coverage in (min_samples) number of samples (i.e. the gene is expressed in (min _samples) number of samples) '''
								open(outfile_parsed,'a').write('%s\t%s\t%s\t%s\t%s\t%s\n' % (key, vcf[key][0][SNP], vcf[key][1][SNP], vcf[key][2][SNP][0], coverage, vcf[key][4][SNP]))										
							else:
								print 'min coverage not reached for all files or the variants are complex!'
								pass
				else:
					pass
		else:
			pass

import sys

try:
	in_yaml = sys.argv[1]
	outfile_parsed= sys.argv[2]
	vcf_files = sys.argv[3]
	bam_file_list = sys.argv[4]
	N = sys.argv[5]
	IC = sys.argv[6]
	min_coverage = sys.argv[7]
	min_samples = sys.argv[8]
except:
	print 'check command line arguments!'



snp_parser(in_yaml,outfile_parsed,vcf_files,bam_file_list,N,IC,min_coverage,min_samples)


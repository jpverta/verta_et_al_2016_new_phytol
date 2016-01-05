#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
Given a list of SNP cordinates (genes and positions), produce a sample-specific table of genotypes in each SNP.
Requires:
samples to be processed: a file with names of sampels one per line, these must match names on .vcf files
output from snp_parser_v2.py 
.vcf files
'''

def position_coverage(sample, chromosome, position, bam_file_list='filtered_bam_files.txt',vcf_type=0):
	import pysam
	coverage = 0
    file_list = open(bam_file_list,'r').read().split('\n')
	for bam in file_list: # for each file in list
	if sample in bam:
		try:
			samfile = pysam.Samfile(bam,'rb') # read BAM file with pysam
			sampileup = samfile.pileup(chromosome, int(position)-1,int(position))
			for pileupcolumn in sampileup:
				if vcf_type == 0: # if SNP positions in vcf file correspond exactly to positions in BAM file
					if pileupcolumn.pos == int(position): # look above concerning cordinates
						coverage = int(pileupcolumn.n)
				elif vcf_type == -1: # if SNP positions in VCF file are shifted -1 relative to BAM file (SAMTOOLS mpileup)
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
    return coverage




def genotypeTable(vcfFile,SNPfile,vcfPath):
	import vcf
	for vcfSample in open(vcfFile,'r').read().split('\n'):
		sample = vcfSample.split('B.vcf')[0]
		outfile = open(sample+'.genot','a')
		outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % ('Gene','Position','Ref', 'Alt', 'Mgg_SNP','Mgg_cov','Embryo_AF','Embryo_Alt', 'Embryo_AO', 'Embryo_RO'))
		for SNPline in open(SNPfile,'r').read().split('\n'): #read snp.out file line by line
			if len(SNPline.split('\t')) > 4:
				mgg_alt = list(SNPline.split('\t')[5].split(', \'')) # the samples for which the SNP exists
				mgg_alt[0] = mgg_alt[0].split('[\'')[1] # remove the brackets fro the first sample on the list
				mgg_alt = [x.split('A')[0] for x in mgg_alt]
				gene = SNPline.split('\t')[0] # the gene with the SNP
				pos = SNPline.split('\t')[1] # position of the SNP
				ref = SNPline.split('\t')[2][2] # reference SNP
				alt = SNPline.split('\t')[3][0] # alternative SNP (observed)
				vcf_reader = vcf.Reader(open(vcfPath+sample+'B.vcf.gz','r')) # read embryo .vcf file corresponding to genotyped sample
				try: # if vcf record exists it means that the embryo sample cannot be homozygous to the reference SNP, there must be another allele, in either heterozygous or homozygous form
					vcfRecord = vcf_reader.fetch(gene,int(pos)) # fetch the corresponding position in the embryo .vcf file
					if sample in mgg_alt:
						outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, pos, ref, alt, alt, position_coverage(sample,gene,pos,bam_file_list='/home/jpverta/data/Megagametophyte/variants/local/filtered_bam_files.txt'),vcfRecord.INFO['AF'], vcfRecord.ALT, vcfRecord.INFO['AO'], vcfRecord.INFO['RO']))		
					else:
						outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, pos, ref, alt, ref,  position_coverage(sample,gene,pos,bam_file_list='/home/jpverta/data/Megagametophyte/variants/local/filtered_bam_files.txt'), vcfRecord.INFO['AF'], vcfRecord.ALT,  vcfRecord.INFO['AO'], vcfRecord.INFO['RO']))			
				except: # if vcf record doesn't exist the embryo sample must be a homozygous for the reference allele
					cov = position_coverage(sample,gene,pos)
					if sample in mgg_alt and cov > 0: # this category cannot exist: the embryo cannot be homozygous ref if the megagametophyte is alt
						outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, pos, ref, alt, alt, position_coverage(sample,gene,pos,bam_file_list='/home/jpverta/data/Megagametophyte/variants/local/filtered_bam_files.txt'), 'Na', 'Na', 'Na', cov))		
					elif sample not in mgg_alt and cov > 0:
						outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, pos, ref, alt, ref, position_coverage(sample,gene,pos,bam_file_list='/home/jpverta/data/Megagametophyte/variants/local/filtered_bam_files.txt'), '[1.0]', 'Na', '[0]', cov))
					else:
						outfile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (gene, pos, ref, alt, ref,  position_coverage(sample,gene,pos,bam_file_list='/home/jpverta/data/Megagametophyte/variants/local/filtered_bam_files.txt'),'Na', 'Na', 'Na', cov))
			else:
				break

import sys 

try:
	vcfFile = sys.argv[1]
	SNPfile = sys.argv[2]
	vcfPath = sys.argv[3]	
except:
	print 'Check command line arguments!'

genotypeTable(vcfFile,SNPfile,vcfPath)


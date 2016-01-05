#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
Script for validating base calls (in .bam files) across samples and SNPs. 
Input:
Text file containing genes and positions to be analyzed, one gene and position per line.
output from vcf_parser.py (.yaml file containing SNP calls)
'''


def position_nucleotides(chromosome, position, bam_file_list):
	import pysam
	import re
	samples =  [re.sub('[\D]','',bam) for bam in open(bam_file_list,'r').read().split('\n') if len(re.sub('[\D]','',bam)) > 1]
    coverage = dict.fromkeys(samples)
    file_list = open(bam_file_list,'r').read().split('\n')
    for sample in samples:
		for bam in file_list: # for each file in list
           	if re.sub('[\D]','',bam) == sample:
				coverage[sample]=[]
				samfile = pysam.Samfile(bam,'rb') # read BAM file with pysam
                sampileup = samfile.pileup(chromosome, int(position))
                for pileupcolumn in sampileup: # look above concerning cordinates
					if pileupcolumn.pos == int(position): # look above concerning cordinates
						for pileupread in pileupcolumn.pileups:
							coverage[sample].append(pileupread.alignment.seq[pileupread.qpos-1])
        return coverage


def genotypeValidator(pos_file,outfile,bam_list):
	import yaml
	import re
	vcf = open(pos_file,'r').read()
	samples =  [re.sub('[\D]','',bam) for bam in open(bam_file_list,'r').read().split('\n') if len(re.sub('[\D]','',bam)) > 1]
	open(outfile,'w').write('Gene'+'\t'+'Pos'+'\t')
	for s in samples:
		open(outfile,'a').write(s+'\t')
		open(outfile,'a').write('\n')
			for line in vcf.split('\n'):
				gene = line.split(' ')[0]
				position = line.split(' ')[1]
				open(outfile,'a').write(str(gene)+'\t'+str(position)+'\t')
				nucleotides = position_nucleotides(str(gene),int(position),bam_list)
				for s in samples:
					if len(nucleotides[s]) > 1:
                  		open(outfile, 'a').write(str(float(sum([N == nucleotides[s][0] for N in nucleotides[s]]))/float(len(nucleotides[s])))+'\t')
                	else:
						open(outfile, 'a').write('NA'+'\t')	
				open(outfile, 'a').write('\n')


import sys
                                        
try:
        pos_file = sys.argv[1]
        outfile_parsed= sys.argv[2]
        bam_file_list = sys.argv[3]
except:
        print 'check command line arguments!'
                                
                                
genotypeValidator(pos_file,outfile_parsed,bam_file_list)


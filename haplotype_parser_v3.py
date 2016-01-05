#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
Script to parse linked SNPs per gene from a file containing all observed SNPs (output from vcf_parser.py -script).

Input:
Yaml -file containing the SNP dictionary (vcf_parser.py output)
Number of samples.
Confidence interval for segregation frequencies IC +- N/2.
Minimum coverage of SNP position required per sample.
Minimum number of samples where the coverage is attained.
File containing paths to .vcf files, one path per line.
File containing paths to .bam files, one path per line.
Output from linkageValidation.py 
 		
'''

def position_type(chromosome, position, vcf_file_list):
	import vcf
	var_type = []
	file_list = open(vcf_file_list,'r').read().split('\n')
	file_list = [file+'.gz' for file in file_list]
	for file in file_list: # for each file in list
		try:	
			reader = vcf.Reader(filename=file)
			for record in reader.fetch(chromosome,position,position+1): # look above concerning cordinates
				var_type.append(record.INFO['TYPE'][0])
		except:
			print 'Warning in position_type: Could not open VCF file %s' % file
			pass
	return set(var_type)



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
				for pileupcolumn in sampileup: # look above concerning cordinates
					if pileupcolumn.pos == position: # look above concerning cordinates
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

def haplotyper(in_dict, gene, positions, vcf_files, samples):
	import re
	import itertools
	N=66
	IC=8
    sample_dict = dict()
	for SNP in  [in_dict[gene][0].index(int(SNP)) for SNP in positions]: # for each index of SNPs in dictionary
		sample_dict[SNP] = [re.sub("[A-z._' ]",'',sample) for sample in in_dict[gene][4][SNP]] # strip alt sample list to contain only numbers
        hap_dict = dict.fromkeys(sample_dict.keys()) # a dictionary that gives the linkage information between SNPs, keys: SNP indexes, values: the other SNP indexes that are in linkage
        SNP_filter = [SNP_entry for SNP_entry in sample_dict if len(sample_dict[SNP_entry]) in range(int(N)/2-int(IC) ,int(N)/2+int(IC)) and sum([qual >= 30 for qual in quality_filter(gene,in_dict[gene][0][SNP_entry],vcf_files)]) >= 0.9*len(sample_dict[SNP_entry])] # filter only the SNP entries that might are observed in a Mendelian$
        for SNP_entry in SNP_filter: # for each SNP in the filtered list
			hap_dict[SNP_entry] = [[],[]]
            for target_entry in SNP_filter: # for each of the other MENDELIAN SNP entries in the list of SNPs in the gene (do not use a genertor at this point!)
					if target_entry is not SNP_entry and len(sample_dict[target_entry]) in range(int(N)/2-int(IC) ,int(N)/2+int(IC)): # excluding self-self comparisons
                       	SNP_samples = sample_dict[SNP_entry]
                        target_samples = sample_dict[target_entry]
                        sample_match = list()
                        sample_match.append([x in target_samples for x in SNP_samples]) # comparison A to B
                        sample_match.append([x in SNP_samples for x in target_samples]) # comparison B to A
                        if sum(sample_match[0]) >= 0.9*len(sample_match[0]) and  sum(sample_match[1]) >= 0.9*len(sample_match[1]):
							''' if the SNPs are observed in same samples 90% of the time, this is a two-way comparison, A in B and B in A '''
							if  not hap_dict[SNP_entry][0]: # if this is the first comaprison for the SNP position
								hap_dict[SNP_entry][0] = [target_entry]
							else: # if not
								hap_dict[SNP_entry][0].append(target_entry)
						elif len(set(SNP_samples+target_samples)) == len(SNP_samples+target_samples) and len(SNP_samples+target_samples) >= 0.9*len(samples): 
							# id each sample is present only once in each list -> the two lists reconstitute the list of all samples
							''' this is a comparison of alternative SNPs relative to the referece: if the sample lists re-constitute the sample list then they are linked 
							--> hence you have to distribute the pair to different lists --> you have to assign the reference status to one set of samples so that 
							when you summarize read counts per linked sets of SNPs the haplotypes have the right counts --> the haplotype will have both alternative AND
							reference counts assigned to it --> test with e.g. WS00749_F07.1'''
							if not hap_dict[SNP_entry][1]: # if this is the first comparison for the SNP position
								hap_dict[SNP_entry][1] = [target_entry]
							else: # if not
								hap_dict[SNP_entry][1].append(target_entry)
	hap_list = [[],[]]
	''' 
	hap_dict: keys are SNP indexes, each value consists of two lists, 1st list gives the other SNPs that are observed in the same samples
	2nd list gives the SNPs that are same as reference in the same samples, but contain a linked SNP in the alternative samples --> the SNPs are still linked:
	INDEX	 0   1   2   3   4
	ALT 1 ---G---A---T---T---G---
	REF   ---G---A---T---A---C---
	ALT2  ---A---C---G---A---C---
	Example, for index 0: 1st list would contain the SNPs 0,1,2
						  2nd list would contain the SNPs 3,4
	--> the lists always give the positions that show the linked, alternative SNPs 
	this gives the two maternal haplotypes:
	ALT1 -G-A-T-T-C-
	ALT2 -A-C-G-A-C- 
	'''
	if len(hap_dict.keys()) == 2 and all([hap_dict[x] for x in hap_dict.keys()]) and all([len(hap_dict[x][0]) == 0 for x in hap_dict.keys()]):
		# this is a special case where two complementary SNPs are observed, i.e. when the alternative SNPs relative to the reference are on the two different alleles
		return [set(hap_dict[hap_dict.keys()[1]][1]),set(hap_dict[hap_dict.keys()[0]][1])]
	for SNP in hap_dict.keys():
		if len(hap_list[0]) > 1 and len(hap_list[1]) > 1: # if there are already SNPs for this gene in hap_list
			if hap_dict[SNP] is None:
				pass 
			elif SNP in hap_list[0]:
				hap_list[0].extend(hap_dict[SNP][0])
				hap_list[1].extend(hap_dict[SNP][1])	
			elif SNP in hap_list[1]:
				hap_list[0].extend(hap_dict[SNP][1])
				hap_list[1].extend(hap_dict[SNP][0])
			elif any([x in hap_list[0] for x in hap_dict[SNP][0]]): # if the SNP itself is not observed in the lists but any of its positive partners are
				hap_list[0].extend(hap_dict[SNP][0])
				hap_list[1].extend(hap_dict[SNP][1])
			elif any([x in hap_list[1] for x in hap_dict[SNP][0]]): 
				hap_list[0].extend(hap_dict[SNP][1])
				hap_list[1].extend(hap_dict[SNP][0])
			else:
				print 'haplotype not resolved for SNP '+ str(SNP)
		elif hap_dict[SNP] and len(hap_dict[SNP][0]) >= 1 and len(hap_dict[SNP][1]) >= 1: # if this is the first SNP for the gene, initiate the lists of linked SNP representing ALT and REF nucleotides
			hap_list[0].extend(hap_dict[SNP][0])
			hap_list[1].extend(hap_dict[SNP][1])
		elif hap_dict[SNP] and len(hap_dict[SNP][0]) >= 1: # 
			hap_list[0].extend(hap_dict[SNP][0])
		elif hap_dict[SNP] and len(hap_dict[SNP][1]) >= 1: # 
			hap_list[0].extend(hap_dict[SNP][0])
		else:
			print 'no linked partners for SNP ' + str(SNP)
	print 'hap_list: '
	print hap_list[0]
	print hap_list[1]
	print set(map(int,hap_list[0]))
	print set(map(int,hap_list[1]))
	if len(hap_list[0]) >= 2 or len(hap_list[1]) >= 2:
        	return [set(map(int,hap_list[0])),set(map(int,hap_list[1]))]
	else:
		return [[],[]]


def get_linkage(gene,linkageFile):
	import re
	for line in open(linkageFile,'r').read().split('\n'):
		if line.split('\t')[0] == gene:
			linkage = dict(zip( [re.sub("[\D]",'',x) for x in line.split('\t')[1].split(',')], [re.sub("[\D]",'',x) for x in line.split('\t')[2].split(',')]))
			return linkage


def snp_parser(in_yaml,outfile_parsed,N,IC,min_coverage,min_samples,vcf_files,bam_file_list,linkageFile='linkageValidation.txt'):
	''' Linkage file determines whether all the SNPs that have been observedfor the gene are in linkage over all samples, indicated as 1 on the third field of each row (transcript). 0 means no linkage
	over samples, which might be caused by reads from another transcript mapping to the same GCAT gene model sequence. These should be evited.
	'''
	print 'running with arguments ' + in_yaml,outfile_parsed,vcf_files,bam_file_list,N,IC,min_coverage,min_samples
	import yaml, re
	open(outfile_parsed,'a')
	bases = []
	vcf = yaml.load(open(in_yaml,'r').read())
	samples =  [re.sub('[\D]','',bam) for bam in open(bam_file_list,'r').read().split('\n') if len(re.sub('[\D]','',bam)) > 1] 
	a_os = dict([(s,0) for s in samples])
	r_os = dict([(s,0) for s in samples])
	print r_os
	print a_os
	for key in vcf.keys(): # for all chromosomes (transcripts)
		linkage = get_linkage(key,linkageFile)
		print key
		print linkage
		if len(vcf[key]) > 1 and linkage and len([x for x in linkage.keys() if linkage[x] == '1']) >= 2 and (sum([x == '1' for x in linkage.values()]) >= 0.9*len(linkage.values())): # if it has SNPs and 90% the SNPs are under linkage
			haplo = haplotyper(vcf, key, [x for x in linkage.keys() if linkage[x] == '1'], vcf_files, samples) # haplotype parsing function
			open('haplotype_lists.out','a').write(key+'\t'+str([vcf[key][0][x] for x in haplo[0]])+'\t'+str([vcf[key][0][x] for x in haplo[1]])+'\n')
			alt_samples = dict()
			ref_samples = dict()
			alt_coverage = dict()
			ref_coverage = dict()
			variant_types = dict()
			for SNP_index in haplo[0]: # for each SNP position in haplotype for the 1st positions that share the ALT SNP 
				alt_samples[vcf[key][0][SNP_index]] = [re.sub('[\D]','',sample) for sample in vcf[key][4][SNP_index]]
				ref_samples[vcf[key][0][SNP_index]] = [s for s in samples if s not in alt_samples[vcf[key][0][SNP_index]]]
				alt_coverage[vcf[key][0][SNP_index]] = position_coverage(key, vcf[key][0][SNP_index], alt_samples[vcf[key][0][SNP_index]], bam_file_list) # record the coverage of the position in all alt samples 
				ref_coverage[vcf[key][0][SNP_index]] = position_coverage(key, vcf[key][0][SNP_index], ref_samples[vcf[key][0][SNP_index]], bam_file_list)  # coverage all ref samples 	
				for s in alt_samples[vcf[key][0][SNP_index]]:
					a_os[s] += alt_coverage[vcf[key][0][SNP_index]][s]
				for s in ref_samples[vcf[key][0][SNP_index]]:
					r_os[s] += ref_coverage[vcf[key][0][SNP_index]][s]
				assert all([s in alt_coverage.keys() for s in alt_samples])
				assert all([s in ref_coverage.keys() for s in ref_samples])
			for SNP_index in haplo[1]: 
				'''
				for each SNP position in haplotype for the 2st positions that share the ALT SNP 
				the samples in the ALT field of the vcf dictionary have to be counted as the the reference to have the coverage of the two maternal haplotypes
				the REF SNP is linked to the ALT SNPs of the positions in the 1st list of "haplo" --> see haplotyper function   
				'''
				ref_samples[vcf[key][0][SNP_index]] = [re.sub('[\D]','',sample) for sample in vcf[key][4][SNP_index]] # these samples actually have the ALT nucleotide relative to the reference! 
				alt_samples[vcf[key][0][SNP_index]] = [s for s in samples if s not in ref_samples[vcf[key][0][SNP_index]]]
				ref_coverage[vcf[key][0][SNP_index]] = position_coverage(key, vcf[key][0][SNP_index], ref_samples[vcf[key][0][SNP_index]], bam_file_list) # record the coverage of the position in all alt samples 
				alt_coverage[vcf[key][0][SNP_index]] = position_coverage(key, vcf[key][0][SNP_index], alt_samples[vcf[key][0][SNP_index]], bam_file_list)  # coverage all ref samples 	
				for s in ref_samples[vcf[key][0][SNP_index]]:
					a_os[s] += ref_coverage[vcf[key][0][SNP_index]][s] # the "REF" reads cover actually ALT nucleotides, and need to be counted to the alternative observations...
				for s in alt_samples[vcf[key][0][SNP_index]]:
					r_os[s] += alt_coverage[vcf[key][0][SNP_index]][s] # and vice versa!				
				assert all([s in alt_coverage.keys() for s in alt_samples])
				assert all([s in ref_coverage.keys() for s in ref_samples])
			if len(haplo[0]) > 0 or len(haplo[1]) > 0 and  [sum([x >= int(min_coverage) for x in alt_coverage[SNP_index].values()]) for SNP_index in alt_coverage.keys()] +  [sum([x >= int(min_coverage) for x in ref_coverage[SNP_index].values()]) for SNP_index in ref_coverage.keys()] >= int(min_samples): 
				open(outfile_parsed,'a').write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (key, [len(vcf[key][4][SNP]) for SNP in haplo[0] | haplo[1]], [vcf[key][0][SNP] for SNP in haplo[0] | haplo[1]], [vcf[key][2][SNP] for SNP in haplo[0] | haplo[1]], [alt_coverage[vcf[key][0][SNP]] for SNP in haplo[0] | haplo[1]], [ref_coverage[vcf[key][0][SNP]] for SNP in haplo[0] | haplo[1]], alt_samples, ref_samples))
				print 'RUN COMPLETED!!!'										
			else:
				pass
		else:
			pass
	for s in samples:
		open('alt_observations.haplotypes','a').write(str(s) + ' ' + str(a_os[s]) + '\n')
		open('ref_observations.haplotypes','a').write(str(s) + ' ' + str(r_os[s]) + '\n')



import sys

try:
	in_yaml = sys.argv[1]
	outfile_parsed= sys.argv[2]
	N = sys.argv[3]
	IC = sys.argv[4]
	min_coverage = sys.argv[5]
	min_samples = sys.argv[6]
	vcf_files = sys.argv[7]
	bam_file_list = sys.argv[8]
except:
	print 'check command line arguments!'


snp_parser(in_yaml,outfile_parsed,N,IC,min_coverage,min_samples,vcf_files,bam_file_list)


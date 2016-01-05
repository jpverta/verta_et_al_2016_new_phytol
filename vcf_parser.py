#!/home/jpverta/bin/python2.7

''' 
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014

SNP DICTIONARY "vcf_parser":

	for all .flt.vcf files in a directory
	pick a set of SNPs that are found segregating 1:1
	store data in a dictionary:
	{chromosomes : positions;SNPs;observations (in samples)}
	for example:
	('GQ03326_A22.1': ([341,[A,T],[14,6]], [566,[G,[C,T]],[10,10]])} multiple SNPs within a chromosome separated by a comma
	
'''

def get_chromosomes(infile):
	with open(infile,'r') as file:
		for line in file:
			yield line.split('\n')[0], ()
			

def vcf_parser(input_files,min_depth,filter_complex,chromosome_file,verbose,yaml_file): #input: text file with a line for each files to analyse, chromosomes, minimum depth for SNP call
	import vcf
	import yaml
	alphabet = {'[A, T]':'W','[T, A]':'W','[C, G]':'S','[G, C]':'S','[A, C]':'M','[C, A]':'M',\
	'[G, T]':'K','[T, G]':'K','[A, G]':'R','[G, A]':'R','[C, T]':'Y','[T, C]':'Y'}
	chromosomes = dict(get_chromosomes(chromosome_file)) #create a dictionary with keys as chromosome names
	file_list = open(input_files,'r').read().split('\n')# list of files to analyse
	for file in file_list: # for each file in list
		try: 
			reader = vcf.Reader(open(file,'r')) #vcf file reader
			for record in reader:	# for each SNP call in vcf file
				if record.INFO['DP'] >= float(min_depth): # if SNP is observed over minimum depth 
					if len(record.ALT) == 1 and record.INFO['TYPE'] == ['snp']: #if the alternative SNP is NOT POLYMORPHIC and a simple SNP (not a series of nucleotides)'''
						if len(chromosomes[record.CHROM]) == 0:	#if a SNP instance doesnt's exist for the cromosome
							chromosomes[record.CHROM] = ([record.POS],[[record.REF]],[[str(record.ALT[0])]],[[1]],[[file]]) # create lists for position, ref SNP, alt SNP and observation
						elif record.POS not in chromosomes[record.CHROM][0]: #if SNP position for a chromosome doesn't exist
							chromosomes[record.CHROM][0].append(record.POS)
							chromosomes[record.CHROM][1].append([record.REF])
							chromosomes[record.CHROM][2].append([str(record.ALT[0])])
							chromosomes[record.CHROM][3].append([1])
							chromosomes[record.CHROM][4].append([file])
						elif record.POS in chromosomes[record.CHROM][0]: #if position exists
							posind = chromosomes[record.CHROM][0].index(record.POS) # index of the position within the values
							if str(record.ALT[0]) not in str(chromosomes[record.CHROM][2][posind]): # if the alternative SNP doesn't exist for the position
								chromosomes[record.CHROM][2][posind].append(str(record.ALT[0])) # append the alternative SNP to ALT list under the position
								chromosomes[record.CHROM][3][posind].append(1) # record the number it has been observed (in samples) - last in list
								chromosomes[record.CHROM][4][posind].append(file)
							elif str(record.ALT[0]) in str(chromosomes[record.CHROM][2][posind]): # if the alternative SNP exists for the position
								ind=chromosomes[record.CHROM][2][posind].index(str(record.ALT[0]))# record it's index
								chromosomes[record.CHROM][3][posind][ind]+=1 # add one to the number it has been observed (in samples)
								chromosomes[record.CHROM][4][posind].append(file)
							else:
								print 'Checkpoint 1'
								break
						else:
							print 'Checkpoint 2'
							break
					elif len(record.ALT) == 2: # if the alternative SNP is POLYMORPHIC -  exists only in vcf files produced with SAMTOOLS'''
						record.ALT = [alphabet[str(record.ALT)]]
						if len(chromosomes[record.CHROM]) == 0:	#if a SNP instance doesnt's exist for the cromosome
							chromosomes[record.CHROM] = ([record.POS],[record.REF],[record.ALT],[[1]],[[file]])
						elif record.POS not in chromosomes[record.CHROM][0]: #if SNP position for a chromosome doesn't exist	
							chromosomes[record.CHROM][0].append(record.POS)
							chromosomes[record.CHROM][1].append(record.REF)
							chromosomes[record.CHROM][2].append(record.ALT)
							chromosomes[record.CHROM][3].append([1])
							chromosomes[record.CHROM][4].append([file])
						elif record.POS in chromosomes[record.CHROM][0]: #if position exists
							posind = chromosomes[record.CHROM][0].index(record.POS) # index of the position within the values
							if str(record.ALT[0]) not in str(chromosomes[record.CHROM][2][posind]): # if the alternative SNP doesn't exist for the position
								chromosomes[record.CHROM][2][posind].append(record.ALT[0]) # append the alternative SNP to ALT list under the position
								chromosomes[record.CHROM][3][posind].append(1) # record the number it has been observed (in samples) - last in list
								chromosomes[record.CHROM][4][posind].append(file)
							elif str(record.ALT[0]) in str(chromosomes[record.CHROM][2][posind]): # if the alternative SNP exists for the position
								ind=chromosomes[record.CHROM][2][posind].index(record.ALT[0])# record it's index
								chromosomes[record.CHROM][3][posind][ind]+=1 # add one to the number it has been observed (in samples)
								chromosomes[record.CHROM][4][posind].append(file)
							else:
								print 'Checkpoint 3'
								break
						else:
							print 'Checkpoint 4'
							break 
					elif verbose == True:
						print 'multiallelic SNP or complex variant in %s' % record
						pass
					else:
						pass
				else:
					pass
		except:
			print record
			print 'Could not process file %s' % file
			pass
	if not yaml_file:
		return chromosomes
	else:
		open(yaml_file,'w').write(yaml.dump(chromosomes))


import sys

try:
	input_files = sys.argv[1]
	min_depth= sys.argv[2]
	filter_complex = sys.argv[3]
	chromosome_file = sys.argv[4]
	verbose = sys.argv[5]
	yaml_file= sys.argv[6]
except:
	print "check commandline arguments!"

vcf_parser(input_files,min_depth,filter_complex,chromosome_file,verbose,yaml_file)

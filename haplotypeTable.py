#!/home/jpverta/bin/python2.7

'''
J-P Verta, Universite Laval, Quebec, Canada, 2013-2014
Given a list of haplotype cordinates (genes and positions), produce a sample-specific table of genotypes in each haplotype (set of linked SNPs within a single gene).
Requires:
samples to be processed: a file with names of sampels one per line, these must match names on .vcf files
output from haplotype_parser_v3.py (haplotype_list.out -file)
output from genotypeTable.py
.vcf files
'''

def haplotyper(infile):
        import ast
        haplofile = open(infile,'r').read()
        haplo = dict.fromkeys([l.split('\t')[0] for l in haplofile.split('\n') if len(l) > 1])
        for line in haplofile.split('\n'):
                if len(line) > 2:
                        haplo[line.split('\t')[0]] = [[],[]]
                        haplo[line.split('\t')[0]][0] = ast.literal_eval(line.split('\t')[1])
                        haplo[line.split('\t')[0]][1] = ast.literal_eval(line.split('\t')[2])
        return haplo



def haplotype_dict(sample,haplolist):
	snpfile = open(sample+'.genot','r').read()
	haplo_groups = haplotyper(haplolist)
	#Fecth all SNP lines that match the haplotype key (gene)
	AF = dict()
	positions = dict()
	hap_A = dict()
	hap_B = dict()
	for gene in haplo_groups.keys():
		SNPlines = []
		if haplo_groups[gene] and len(haplo_groups[gene][0]) > 0 or len(haplo_groups[gene][1]) > 0: # if the gene has a haplotype 
			AF[gene] = []
			positions[gene] = []
			hap_A[gene] = []
			hap_B[gene] = []
			for SNPline in snpfile.split('\n'):
				if SNPline.split('\t')[0] == gene:
					SNPlines.append(SNPline)
			for SNP in SNPlines:
				pos = SNP.split('\t')[1]
				positions[gene].append(pos)
				AF[gene].append(SNP.split('\t')[6])
				if int(pos) in haplo_groups[gene][0]:
					hap_A[gene].append(SNP.split('\t')[8]) # record the 'AO' field into the hap_A list
					hap_B[gene].append(SNP.split('\t')[9]) # record the 'RO' field into the hap_B list
				elif int(pos) in haplo_groups[gene][1]:
					hap_A[gene].append(SNP.split('\t')[9]) # vice versa
					hap_B[gene].append(SNP.split('\t')[8])
				else:
					print 'cannot find SNP position in haplolist: ' + sample + ' ' + SNP  
					pass
	return [positions, AF, hap_A, hap_B]


def run(vcfFiles,haplolist):
	import os, re
	for vcfSample in open(vcfFiles,'r').read().split('\n'):
		sample = vcfSample.split('B.vcf')[0]
		if sample+'.haplot' not in os.listdir('.'):
			outfile = open(sample+'.haplot','a')
            outfile.write('%s\t%s\t%s\t%s\t%s\n' % ('Gene','Positions','AF','Embryo_A_obs','Embryo_B_obs'))
            haplotypes = haplotype_dict(sample,haplolist)
			print haplotypes
			print haplotypes[0].keys()
			for gene in haplotypes[0].keys():
				outfile.write('%s\t%s\t%s\t%s\t%s\n' % (gene, [pos for pos in haplotypes[0][gene]], [af for af in haplotypes[1][gene]], [ao for ao in haplotypes[2][gene]], [ro for ro in haplotypes[3][gene]]))
		print 'Run completed for sample ' + str(sample)



import sys

try:
        vcfFiles = sys.argv[1]
        haplolist = sys.argv[2]
except:
        print 'Check command line arguments!'

run(vcfFiles, haplolist)

	

# R version 2.15.2

'''
Loop throught all .genot files in a directory and summarize the genotype observations to a table. Output genotype frequencies for SNP position that are not covered by a haplotype.
Requires:
output from haplotypeTable.py (.genot files)
output from genotypeTable.py (.haplot files)
'''

genoPath = # PATH TO A SINGLE .genot FILE
haploPath = # PATH TO A SINGLE .haplot FILE
outPath = # PATH TO OUTPUT FILE, FOR EXAMPLE: singlesTable.txt

library(stringr)

aseTable = read.table(genoPath,header=T,sep='\t')
aseTable = aseTable[,1:2]
aseTable[,3:6] = 0
rownames(aseTable) = paste(aseTable[,1],aseTable[,2],sep='_')

haploTable = read.table(haploPath,header=T,sep='\t')
rownames(haploTable)=haploTable[,1] 

NAtoZero = function(x){
	x[is.na(x)] = 0
	return(x)
}

for (i in list.files()){
	if (grepl('.genot',i)==T){
		genot = read.table(i,header=T,sep='\t')
		rownames(genot) = paste(genot[,1],genot[,2],sep='_')
		singles = genot[genot[,1] %in% haploTable[,1] == F,]
		singleGenes = as.character(singles[,1])
		singleMatch1 = aseTable[,1] %in% singleGenes
		#all samples:
		aseTable[singleMatch1,6] = as.numeric(aseTable[singleMatch1,6]) + 1
		#heterozygous samples:
		singles = singles[singles[,7] == '[0.5]',]
		singleGenes = as.character(singles[,1])
		singleMatch = aseTable[,1] %in% singleGenes
		aseTable[singleMatch,3] = as.numeric(aseTable[singleMatch,3]) + 1
		aseTable[singleMatch,4] = as.numeric(aseTable[singleMatch,4]) + NAtoZero(as.numeric(str_replace_all(singles[singles[,1] %in% aseTable[singleMatch,1],9],"[[:punct:]]","")))
		aseTable[singleMatch,5] = as.numeric(aseTable[singleMatch,5]) + NAtoZero(as.numeric(as.character(singles[singles[,1] %in% aseTable[singleMatch,1],10])))
	}
}

colnames(aseTable)[3:6] = c('Heteroz_obs','Embryo_AO','Embryo_BO','All_obs')

write.table(aseTable[singleMatch1,],outPath)

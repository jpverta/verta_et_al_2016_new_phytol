# R version 2.15.2

'''
Loop through all .haplot files in a directory and summarize the genotype observations to a table.
Requires:
output from haplotypeTable.py (.haplot files)
'''

library(stringr)

haploPath = # PATH TO A SINGLE .haplot FILE
outPath = # PATH TO OUTPUT FILE, FOR EXAMPLE: haploTable.txt

haploTable = read.table(haploPath,header=T,sep='\t')
haploTable = haploTable[,1:2]
head(haploTable)
haploTable[,3:5] = 0

for (i in list.files()){
	if (grepl('.haplot',i)==T & grepl('.haplotypes',i)==F){
		haplot = read.table(i,header=T,sep='\t')
		for (y in 1:nrow(haplot)){
			if (is.na(haplot[y,3]) == FALSE & (length(which(grepl('0.5', unlist(strsplit(as.character(haplot[y,3]),','))))) > length(unlist(strsplit(as.character(haplot[y,3]),',')))/2)){ 
				# consider the sample as heterozygous if heterozygous SNPs are more frequent in the haplotype than homozygous SNPs 
				haploTable[y,3] = haploTable[y,3]+1
				haploTable[y,4] = haploTable[y,4]+sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplot[y,4]),"[[:punct:]]",""),' ')))))
				haploTable[y,5] = haploTable[y,5]+sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplot[y,5]),"[[:punct:]]",""),' ')))))			}
		}
	}
}

colnames(haploTable)[3:5] = c('Heteroz_obs','Embryo_AO','Embryo_BO')
haploTable[,6] = haploTable[,4]/haploTable[,3]
haploTable[,7] = haploTable[,5]/haploTable[,3]
colnames(haploTable)[6:7] = c('Embryo_AO_mean','Embryo_BO_mean')
write.table(haploTable,haploPath,row.names=F)

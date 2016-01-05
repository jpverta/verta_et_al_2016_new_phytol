# R version 2.15.2

'''
Summarize observed genotype frequencies across .genot files.
Requires:
output from genotypeTable.py (.genot files)
'''

outPath = # PATH TO OUTFILE, FOR EXAMPLE: genoTable.txt

genoPath = # PATH TO A SINGLE .genot FILE, THE SNP IN THIS FILE WILL BE SUMMARIZED

library(stringr)

genoTable = read.table(genoPath,header=T,sep='\t')
genoTable = genoTable[,1:4]
genoTable[,5:7] = 0

for (i in list.files()){
	if (grepl('.genot',i)==T){
		genot = read.table(i,header=T,sep='\t')
		for (y in 1:nrow(genot)){
			if (is.na(genot[y,7]) == FALSE){
				if (is.na(genot[y,7]) == FALSE & genot[y,7] == "[0.5]"){ 
				#the embryo is heterozygous if AF == 0.5
					genoTable[y,6] = genoTable[y,6]+1
				}
				if (genot[y,7] == "[1.0]" & str_replace_all(genot[y,8],"[[:punct:]]","") == "Na" & genot[y,5] %in% genot[y,3] ==F){ 
				#the embryo is heterozygous if no alternative SNP is observed but the mgg genotype not ref 
				#--> at least ONE allele should be the same as in the megagametophyte
					genoTable[y,6] = genoTable[y,6]+1
				}
				if (genot[y,7] == "[1.0]" & str_replace_all(genot[y,8],"[[:punct:]]","") == genot[y,4] & genot[y,5] == genot[y,3]){ 
				#the embryo is heterozygous if only the alternative SNP is observed but the mgg genotype is ref 
				#--> at least ONE allele should be the same as in the megagametophyte
					genoTable[y,6] = genoTable[y,6]+1
				}
				if (genot[y,7] == "[1.0]" & str_replace_all(genot[y,8],"[[:punct:]]","") == "Na" & genot[y,5] == genot[y,3] & genot[y,10] > 0){
				# embryo is homozygous ref, ref observed in mgg
					genoTable[y,5] = genoTable[y,5]+1
				}
				if (genot[y,7] == "[1.0]" & str_replace_all(genot[y,8],"[[:punct:]]","") == genot[y,4] & genot[y,5] == genot[y,4]){
				# embryo is homozygous alt, alt observed in mgg
					genoTable[y,7] = genoTable[y,7]+1
				}
			}
		}
	}
}

colnames(genoTable)[5:7] = c("Homoz_ref_obs",'Heteroz_obs','Homoz_alt_obs')

write.table(genoTable,outPath)


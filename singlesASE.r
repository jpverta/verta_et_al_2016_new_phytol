# R version 1.15.2

'''
ASE test for SNP positions that are covered by single SNPs and not by haplotypes.
Requires:
output from genotypeTable.py (.genot files)
output from haplotypeTable.py (.haplot files)
output from singlesTable.r
'''

library(stringr)
library(edgeR)

singlesPath = # PATH TO OUTPUT FILE OF singlesTable.r
haploPath = # PATH TO A SINGLE .haplot FILE
outPath = # PATH TO A OUTFILE

aseTable = read.table(singlesPath,header=T)
aseTable = aseTable[,1:5]

haploTable = read.table(haploPath,header=T,sep='\t')
rownames(haploTable)=haploTable[,1]

NAtoZero = function(x){
        x[is.na(x)] = 0
        return(x)
}

AO=list()
RO=list()
HZsamples=list()
A_os=list()
R_os=list()
samples={}
for (i in list.files()){
	if (grepl('.genot',i)==T){
		genot = read.table(i,header=T,sep='\t')
		samples = append(samples,unlist(strsplit(i,'\\.'))[1])
		rownames(genot) = paste(genot[,1],genot[,2],sep='_')
		singles = genot[genot[,1] %in% haploTable[,1] == F,]
        singles = singles[singles[,7] == '[0.5]',]
        singleGenes = as.character(singles[,1])
        singleMatch = aseTable[,1] %in% singleGenes
		for (y in 1:nrow(aseTable)){
			if (as.character(aseTable[y,1]) %in% singleGenes){ 
				if (is.na(aseTable[y,3]) == F & aseTable[y,3]>=25 & aseTable[y,3]<=41){ 
					y = rownames(aseTable)[y]
					HZsamples[[y]] = append(HZsamples[[y]],unlist(strsplit(i,'\\.'))[1])
					AO[[y]] = append(AO[[y]],NAtoZero(as.numeric(str_replace_all(singles[y,9],"[[:punct:]]",""))))
					RO[[y]] = append(RO[[y]],NAtoZero(as.numeric(as.character(singles[y,10]))))
					A_os[[unlist(strsplit(i,'\\.'))[1]]] = append(A_os[[unlist(strsplit(i,'\\.'))[1]]],NAtoZero(as.numeric(str_replace_all(singles[y,9],"[[:punct:]]",""))))
					R_os[[unlist(strsplit(i,'\\.'))[1]]] = append(R_os[[unlist(strsplit(i,'\\.'))[1]]],NAtoZero(as.numeric(as.character(singles[y,10]))))
				}
			}
		}
	}
}

aseTable[,6:10] = NA

for (i in names(AO)){
	if (aseTable[i,3] <= 41 & aseTable[i,3] >= 25){
	print(i)
	A = unlist(AO[i])
	R = unlist(RO[i])	
	tcounts = c(A,R)
	tcounts = t(as.matrix(tcounts))
	colnames(tcounts) = 1:ncol(tcounts)
	rownames(tcounts) = i	
	tgroups = {}
    tgroups[1:length(A)] = "G1"
    tgroups[(length(A)+1):(length(A)+length(R))] = "G2"
    tdesign=model.matrix(~tgroups)
    tcounts = DGEList(counts=tcounts,group=tgroups)
    #option 1: offset calculated as the independent sum of alternative and the sum of reference observations, the read counts of alt and ref alleles are offsetted with these numbers, i.e. the alt read count of a specific gene is offsetted against all the alt reads - saved as offset1
	os = log(unlist(c(lapply(A_os[which(names(A_os) %in% unlist(HZsamples[i]))],sum) , lapply(R_os[which(names(R_os) %in% unlist(HZsamples[i]))],sum) )))
		if (length(tcounts$counts[tcounts$counts > 0]) > 0.9*length(tcounts$counts)){ #to test genes that have counts in more than 90% of the samples
		tcounts = estimateGLMCommonDisp(tcounts,tdesign,offset=os)
		tcounts = estimateGLMTagwiseDisp(tcounts,tdesign,offset=os)
		tfit = glmFit(tcounts,tdesign,offset=os)
		tlrt = glmLRT(tfit)
		# NULL MODEL
        rdesign=tdesign
        rdesign[,2]=1
        rcounts = DGEList(counts=tcounts,group=factor(rep(1,times=ncol(tcounts))))
        rcounts = estimateGLMCommonDisp(rcounts,design=NULL,offset=os)
        rcounts = estimateGLMTagwiseDisp(rcounts,design=NULL,offset=os)
        rfit = glmFit(rcounts,design=NULL,offset=os)
        rfit$'design' = rdesign
        rownames(rfit$'design') = colnames(rfit$'counts')
        rlrt = glmLRT(rfit,coef=1)
        rdisp = rlrt$'dispersion'
        R2 = 1-exp(1)^((tlrt$'table'[i,3])/(-66))
		aseTable[i,6] = tlrt$'table'$'logFC'
		aseTable[i,7] = tlrt$'table'$"PValue" 
		aseTable[i,8] = tlrt$'dispersion'[1]
		aseTable[i,9] = rdisp
		aseTable[i,10] = R2
	}
}
}

colnames(aseTable)[6] = 'logFC'
colnames(aseTable)[7] = "exactTest_pvalue"
colnames(aseTable)[8] = 'model_dispersion'
colnames(aseTable)[9] = 'null_dispersion'
colnames(aseTable)[10] = 'pseudo-R2'

write.table(aseTable,outPath)

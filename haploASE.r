# R version 2.15.2

'''
LRT test for ASE based on haplotype-wise counts in heterozygous genotypes.
Requires:
output from haplotypeHeterozygoteFrequencies.r (genes to be tested and their heterozygote frequencies)
output from haplotypeTable.py (.haplot files)
'''

haploPath = # PATH TO HAPLOTYPE TABLE (output from haplotypeHeterozygoteFrequencies.r)
outPath = # PATH TO OUTFILE

library(stringr)
library(edgeR)

haploTable = read.table(haploPath,header=T)
rownames(haploTable) = haploTable[,1]



AO=list()
RO=list()
HZsamples=list()
A_os=list()
R_os=list()
samples={}
positions=list()
for (i in list.files()){
	if (grepl('.haplot',i)==T & grepl('snp',i) == F & grepl('test',i) == FALSE){
		print(i)
		haplo = read.table(i,header=T,sep='\t')
		rownames(haplo) = haplo[,1]
		samples = append(samples,unlist(strsplit(i,'\\.'))[1])
		for (y in rownames(haplo)){
			if (length(which(grepl('0.5', unlist(strsplit(as.character(haplo[y,3]),','))))) > length(unlist(strsplit(as.character(haplo[y,3]),',')))/2){
				pos = as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,2]),"[[:punct:]]",""),' ')))
				if (length(pos) > 1){
					posSet=c()
					for (h in 1:length(pos)){
						if (length(posSet)==0){
							posSet = append(posSet,pos[h])
						}
						else{
							dist = c()
							for (z in 1:length(posSet)){
								dist = append(dist,abs(as.numeric(pos[h]) - as.numeric(posSet[z])))
							}
							if (all(dist>100)){
								posSet = append(posSet,pos[h])
							}
						}
					}
					posSet = posSet[posSet %in% 0 == F] 
					positions[[y]] = append(positions[[y]],posSet)
					positions[[y]] = unique(positions[[y]])
					posInd = which(pos %in% posSet)
					AO[[y]] = append(AO[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,4]),"[[:punct:]]",""),' ')))[posInd])))
					RO[[y]] = append(RO[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,5]),"[[:punct:]]",""),' ')))[posInd])))
					HZsamples[[y]] = append(HZsamples[[y]],unlist(strsplit(i,'\\.'))[1])
					if (haploTable[y,'Heteroz_obs']>=25 & haploTable[y,'Heteroz_obs']<=41){
A_os[[unlist(strsplit(i,'\\.'))[1]]] = append(A_os[[unlist(strsplit(i,'\\.'))[1]]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,4]),"[[:punct:]]",""),' ')))[posInd])))
R_os[[unlist(strsplit(i,'\\.'))[1]]] = append(R_os[[unlist(strsplit(i,'\\.'))[1]]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,5]),"[[:punct:]]",""),' ')))[posInd])))
					}
				}
				else {
					AO[[y]] = append(AO[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,4]),"[[:punct:]]",""),' '))))))
					RO[[y]] = append(RO[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,5]),"[[:punct:]]",""),' '))))))
					HZsamples[[y]] = append(HZsamples[[y]],unlist(strsplit(i,'\\.'))[1])
					if (haploTable[y,'Heteroz_obs']>=25 & haploTable[y,'Heteroz_obs']<=41){
A_os[[unlist(strsplit(i,'\\.'))[1]]] = append(A_os[[unlist(strsplit(i,'\\.'))[1]]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,4]),"[[:punct:]]",""),' '))))))
R_os[[unlist(strsplit(i,'\\.'))[1]]] = append(R_os[[unlist(strsplit(i,'\\.'))[1]]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[y,5]),"[[:punct:]]",""),' '))))))
					}
				}	
			}		
		}
	}
}

# ASE calculation

tryCatch.W.E <- function(expr)
{
    W <- NULL
    w.handler <- function(w){ # warning handler
        W <<- w
        invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),
                                     warning = w.handler),
         warning = W)
}


lrt_test = function(i){
	A = unlist(AO[i])
	R = unlist(RO[i])	
	tcounts = c(A,R)
	tcounts = t(as.matrix(tcounts))
	colnames(tcounts) = 1:ncol(tcounts)
	rownames(tcounts) = i	
	tgroups = {}
	tgroups[1:length(A)] = "Alt"
	tgroups[(length(A)+1):(length(A)+length(R))] = "Ref"
	tdesign=model.matrix(~tgroups)
	tcounts = DGEList(counts=tcounts,group=tgroups)
	os = log(unlist(c(lapply(A_os[which(names(A_os) %in% unlist(HZsamples[i]))],sum) , lapply(R_os[which(names(R_os) %in% unlist(HZsamples[i]))],sum) )))
	tcounts = estimateGLMCommonDisp(tcounts,tdesign,offset=os)
	tcounts = estimateGLMTagwiseDisp(tcounts,tdesign,offset=os)
	tfit = glmFit(tcounts,tdesign,offset=os)
	tlrt = glmLRT(tfit)
	print(tlrt)
	tdisp = tlrt$'dispersion'
	# NULL MODEL
        rdesign=tdesign
        rdesign[,2]=1
        rcounts = DGEList(counts=tcounts,group=factor(rep(1,times=ncol(tcounts))))
        rcounts = estimateGLMCommonDisp(rcounts,design=NULL,offset=os)
        rcounts = estimateGLMTrendedDisp(rcounts,design=NULL,offset=os)
        rcounts = estimateGLMTagwiseDisp(rcounts,design=NULL,offset=os)
        rfit = glmFit(rcounts,design=NULL,offset=os)
		rfit$'design' = rdesign
        rownames(rfit$'design') = colnames(rfit$'counts')
        rlrt = glmLRT(rfit,coef=1)
        rdisp = rlrt$'dispersion'
        R2 = 1-exp(1)^((tlrt$'table'[i,3])/(-66))
	return(c(tlrt$'table',tdisp,rdisp,R2))

}


haploTable[,8] = NA
haploTable[,9] = NA
haploTable[,10] = NA
haploTable[,11] = NA
haploTable[,12] = NA

results=list()
errors = list()
for (i in rownames(haploTable)){
	if (haploTable[i,'Heteroz_obs']>=25 & haploTable[i,'Heteroz_obs']<=41){
		results[i] = list(tryCatch.W.E(lrt_test(i))$value)
		errors[i] = list(tryCatch.W.E(lrt_test(i))$warning)
}
}

for (i in names(results)){
	if (length(results[i][[1]]) > 2){
		haploTable[i,2] = paste(as.character(unlist(positions[i])),collapse=':')
		haploTable[i,8] = results[i][[1]][[1]]
		haploTable[i,9] = results[i][[1]][[4]]
		haploTable[i,10] = results[i][[1]][[5]]
		haploTable[i,11] = results[i][[1]][[6]]
		haploTable[i,12] = results[i][[1]][[7]]	
		
	}
}


colnames(haploTable)[8] = 'logFC'
colnames(haploTable)[9] = "LRT_pvalue"
colnames(haploTable)[10] = "model_dispersion"
colnames(haploTable)[11] = "null_dispersion"
colnames(haploTable)[12] = "pseudo-R2"

write.table(haploTable,outPath)

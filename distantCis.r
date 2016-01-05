'''
J-P Verta 2013-2014

Test for distant effects acting in cis.

Requires:
Output from haplotypeHeterozygoteFrequencies.py
Output from haplotype_parser_v3.py
Significative associations from distantLRT.r
'''

haploPath = # output from haplotypeHeterozygoteFrequencies.py 

variantPath = # output from haplotype_parser_v3.py 

distPath = # output from distantLRTsummary.r 

outPath = # path to output file 

library(stringr)
library(edgeR)
haplo = read.table(haploPath,header=T)
snp = read.table(variantPath,sep='\t')
rownames(snp) = snp[,1]
dist = read.table(distPath,header=T,stringsAsFactors=F)
dist=dist[is.na(dist[,7])==F,]
effects = dist[,2]
genes=dist[,1]

u_genes = unique(genes)
focal = u_genes[u_genes %in% snp[,1]]
haplo_M = haplo[haplo[,'Heteroz_obs'] < 42 & haplo[,'Heteroz_obs'] > 24,]
shared = focal[focal %in% haplo_M[,1]]

# samples that are heterozygous for both focal and distant loci

AO=list()
RO=list()
A_os=list()
R_os=list()
HZsamples=list()
samples={}
for (i in list.files()){
        if (grepl('.haplot',i)==T & grepl('snp',i) == F & grepl('test',i) == FALSE){
                haplo = read.table(i,header=T,sep='\t')
                rownames(haplo) = haplo[,1]
                samples = append(samples,unlist(strsplit(i,'\\.'))[1])
		for (x in 1:length(effects)){
			if (length(which(grepl('0.5', unlist(strsplit(as.character(haplo[effects[x],3]),','))))) > length(unlist(strsplit(as.character(haplo[effects[x],3]),',')))/2 &
			   length(which(grepl('0.5', unlist(strsplit(as.character(haplo[genes[x],3]),','))))) > length(unlist(strsplit(as.character(haplo[genes[x],3]),',')))/2){
                			HZsamples[[rownames(dist)[x]]] = append(HZsamples[[rownames(dist)[x]]],unlist(strsplit(i,'\\.'))[1])
                	}
        	}
	}
}


for (i in names(HZsamples)){
	gene = dist[i,1] #unlist(strsplit(i,':'))[1]
	effect = dist[i,2] #unlist(strsplit(i,':'))[2]
	for (y in unlist(HZsamples[[i]])){
		haplo = read.table(paste(y,'.haplot',sep=''),header=T,sep='\t')
                rownames(haplo) = haplo[,1]
                samples = append(samples,y)
		G_gene=c()
		G_effect=c()
		if (grepl(y,snp[gene,7])){
			G_gene ='A'
		}
		else {
			G_gene = 'B'
		}
 		if (grepl(y,snp[effect,7])){
			G_effect ='A'
		}
		else {
			G_effect = 'B'
		}
		if (G_gene == G_effect){
			AO[[i]] = append(AO[[i]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,4]),"[[:punct:]]",""),' '))))))
                        RO[[i]] = append(RO[[i]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,5]),"[[:punct:]]",""),' '))))))
			A_os[[y]] = append(A_os[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,4]),"[[:punct:]]",""),' '))))))
			R_os[[y]] = append(R_os[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,5]),"[[:punct:]]",""),' '))))))
		}
		else {
			AO[[i]] = append(AO[[i]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,5]),"[[:punct:]]",""),' '))))))
                        RO[[i]] = append(RO[[i]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,4]),"[[:punct:]]",""),' '))))))
			A_os[[y]] = append(A_os[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,5]),"[[:punct:]]",""),' '))))))
			R_os[[y]] = append(R_os[[y]],sum(na.omit(as.numeric(unlist(strsplit(str_replace_all(as.character(haplo[gene,4]),"[[:punct:]]",""),' '))))))
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
        print(i)
        A = unlist(AO[i])
        R = unlist(RO[i])
        print(A)
        print(R)
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


results=list()
errors = list()
for (i in names(HZsamples)){
                results[i] = list(tryCatch.W.E(lrt_test(i))$value)
                errors[i] = list(tryCatch.W.E(lrt_test(i))$warning)
}

resTable=data.frame()
for (i in names(results)){
        if (length(results[i][[1]]) > 2){
                resTable[i,1] = results[i][[1]][[1]]
                resTable[i,2] = results[i][[1]][[4]]
                resTable[i,3] = results[i][[1]][[5]]
                resTable[i,4] = results[i][[1]][[6]]
                resTable[i,5] = results[i][[1]][[7]]
		resTable[i,6] = length(HZsamples[[i]])
        }
}

colnames(resTable)[1] = 'logFC'
colnames(resTable)[2] = "LRT_pvalue"
colnames(resTable)[3] = "model_dispersion"
colnames(resTable)[4] = "null_dispersion"
colnames(resTable)[5] = "pseudo-R2"
colnames(resTable)[6] = "N_double_hetroz"

write.table(resTable,outPath)
		

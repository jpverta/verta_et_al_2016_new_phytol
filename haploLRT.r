
'''
LRT test for differential expression of maternal alleles based on haplotype-wise counts.
'''

haploPath = # PATH TO PARSED HAPLOTYPE COVERAGE -FILE (output from haplotypeCoverageParser.py)

library(edgeR)

haplo=read.table(haploPath,sep='\t',na.strings='None',colClasses = c(rep('character',times=4),rep('numeric',times=132)),skip=1)
haplonames = read.table(haploPath,sep='\t',na.strings='None',nrows=1)
colnames(haplo) = apply(haplonames[1,],2,as.character)
rownames(haplo) = haplo[,1]

alt_counts = read.table('alt_observations.haplotypes',sep='\t',row.names=NULL) 
ref_counts = read.table('ref_observations.haplotypes',sep='\t ',row.names=NULL)

samples = unique(alt_counts[,1])

alt_os =c() 
for (i in as.character(samples)){
	alt_os[i] = sum(alt_counts[alt_counts[,1] %in% i,2])
}

ref_os =c() 
for (i in as.character(samples)){
	ref_os[i] = sum(ref_counts[ref_counts[,1] %in% i,2])
}


for (i in 1:length(alt_os)){
	names(alt_os)[i] = paste(names(alt_os)[i],'_altcov',sep='')
}

for (i in 1:length(ref_os)){
	names(ref_os)[i] = paste(names(ref_os)[i],'_refcov',sep='')
}

alt_os_sums = alt_os+ref_os
ref_os_sums = ref_os+alt_os


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
	tgroups = {}
	altcounts = as.numeric(haplo[i,grepl('alt',colnames(haplo))][haplo[i,grepl('alt',colnames(haplo))] %in% 'NA' == F])
	refcounts = as.numeric(haplo[i,grepl('ref',colnames(haplo))][haplo[i,grepl('ref',colnames(haplo))] %in% 'NA' == F])
	tcounts = c(altcounts,refcounts)
	tcounts = t(as.matrix(tcounts))
	tgroups = c(rep('ALT',times=length(altcounts)),rep('REF',times=length(refcounts)))
	tdesign=model.matrix(~tgroups)
	tcounts = DGEList(counts=tcounts,group=tgroups)
	altnames = names(haplo[i,grepl('alt',colnames(haplo))][haplo[i,grepl('alt',colnames(haplo))] %in% 'NA' == F])
	refnames = names(haplo[i,grepl('ref',colnames(haplo))][haplo[i,grepl('ref',colnames(haplo))] %in% 'NA' == F])
	colnames(tcounts$counts) = c(altnames,refnames)
	a_os = log(alt_os_sums[which(names(alt_os_sums) %in% colnames(tcounts$counts))])
	r_os = log(ref_os_sums[which(names(ref_os_sums) %in% colnames(tcounts$counts))])
	os = c(a_os,r_os)
	tcounts = estimateGLMCommonDisp(tcounts,tdesign,offset = os)
	tcounts = estimateGLMTagwiseDisp(tcounts,tdesign,offset = os)
	tfit = glmFit(tcounts,tdesign,offset = os)
	tlrt = glmLRT(tfit)
	tdisp = tlrt$'dispersion'
	# NULL MODEL
	rdesign=tdesign
	rdesign[,2]=1
	rcounts = DGEList(counts=tcounts,group=factor(rep(1,times=ncol(tcounts$counts))))
	rcounts = estimateGLMCommonDisp(rcounts,design=NULL,offset = os)
	rcounts = estimateGLMTagwiseDisp(rcounts,design=NULL,offset = os)
	rfit = glmFit(rcounts,design=NULL,offset = os)	
	rfit$'design' = rdesign[,]
	rownames(rfit$'design') = colnames(rfit$'counts')	
	rlrt = glmLRT(rfit,coef=1)
	rdisp = rlrt$'dispersion'
	dispRatio = 1-sqrt(tdisp)/sqrt(rdisp)
	R2 = 1-exp(1)^((tlrt$'table'[,3])/(-66))
	return(c(tlrt$'table',tdisp,rdisp,dispRatio,R2))
}

results=list()
for (i in rownames(haplo)){
	print(unlist(lrt_test(i)))
	results[[i]] = unlist(lrt_test(i))
}

resMat = matrix(nrow=length(results),ncol=8)
rownames(resMat) = names(results)
for (i in names(results)){
	if (length(unlist(results[[i]])) == 8){
		resMat[i,1:8] = unlist(results[i])
}
}

colnames(resMat)[1:8]=c('logFC','logCPM','LR','PValue','ModelDispersion','NullDispersion','1-CVRatio','pseudo-R2')	

write.table(resMat,outPath,row.names=TRUE)




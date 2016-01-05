'''
J-P Verta 2013-2014

Local association based on SNP-wise counts. One test for a single SNP.

Requires:
Output from singlesCoverageParser.py
'''

singlePath = # output from singlesCoverageParser.py 

outPath = # path to output file 

library(edgeR)

singles=read.table(singlePath,sep='\t',na.strings='None',colClasses = c(rep('character',times=4),rep('numeric',times=132)),skip=1)
head(singles)
singlesnames = read.table(singlePath,sep='\t',na.strings='None',nrows=1)
singlesnames[1,]
colnames(singles) = apply(singlesnames[1,],2,as.character)
rownames(singles) = singles[,1]

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
	altcounts = as.numeric(singles[i,grepl('alt',colnames(singles))][singles[i,grepl('alt',colnames(singles))] %in% 'NA' == F])
	refcounts = as.numeric(singles[i,grepl('ref',colnames(singles))][singles[i,grepl('ref',colnames(singles))] %in% 'NA' == F])
	tcounts = c(altcounts,refcounts)
	tcounts = t(as.matrix(tcounts))
	tgroups = c(rep('ALT',times=length(altcounts)),rep('REF',times=length(refcounts)))
	tdesign=model.matrix(~tgroups)
	tcounts = DGEList(counts=tcounts,group=tgroups)
	altnames = names(singles[i,grepl('alt',colnames(singles))][singles[i,grepl('alt',colnames(singles))] %in% 'NA' == F])
	refnames = names(singles[i,grepl('ref',colnames(singles))][singles[i,grepl('ref',colnames(singles))] %in% 'NA' == F])
	colnames(tcounts$counts) = c(altnames,refnames)
	a_os =  log(unlist(lapply(apply(apply(singles[,grepl('alt',colnames(singles))],2,as.numeric),2,na.omit),sum)))
	r_os =  log(unlist(lapply(apply(apply(singles[,grepl('ref',colnames(singles))],2,as.numeric),2,na.omit),sum)))
	os = c(a_os,r_os)[which(c(names(a_os),names(r_os)) %in% colnames(tcounts$counts))]
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
for (i in rownames(singles)){
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
resMat	

colnames(resMat)[1:8]=c('logFC','logCPM','LR','PValue','ModelDispersion','NullDispersion','1-CVRatio','pseudo-R2')	

write.table(resMat,outPath,row.names=TRUE)


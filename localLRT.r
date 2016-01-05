'''
LRT test for local effects (differential expression that cosegregates with SNP genotypes) of maternal alleles in haploid megagametophytes.
'''

snpPath = # PATH TO haplotype_parser_v3.py OUTPUT
haploPath = # PATH TO snp_parser_v2.py OUTPUT

'''
Lines (SNPs) from snp file and haplotype file will be parsed to select a non-redundant set of SNP genotypes. Preference will be give for genotypes as defined in haplotype file.
'''

outPath = # PATH TO OUTFILE

countPath = # PATH TO COUNT DATA (output from htseq_count_summary.r), alternatively use the code below to summarize counts across HTSeq output files (.htseq_counts)

'''
counts=data.frame()
for (i in list.files(".")){
if (grepl(".htseq_counts", i)==T){
ht = read.table(i,row.names=1)
counts[1:nrow(ht),i] = ht[,1]
rownames(counts) = rownames(ht)
}
}
'''

library(edgeR)

counts=read.table(countPath)

snp = read.table(snpPath,sep='\t',header=FALSE)
haplo = read.table(haploPath,sep='\t')
rownames(haplo) = haplo[,1]

genes = unique(c(as.character(snp[,1]),as.character(haplo[,1])))

geno = matrix(nrow=length(genes),ncol=1)
rownames(geno) = genes

get_samples = function(x){
	 return(unlist(strsplit(unlist(strsplit(as.character(x),'],'))[2],'\\['))[2])
	}

for (i in rownames(geno)){
	if (i %in% rownames(haplo)){
		geno[i,1] = as.character(get_samples(haplo[i,7])) 
	}
	else {
		if (i %in% snp[,1]){
			geno[i,1] = as.character(snp[,6][snp[,1] %in% i][1])
			
		}	
	}
}


all_counts = DGEList(counts=counts)
all_counts = calcNormFactors(all_counts,method='TMM')

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
	tcounts = counts[i,]
	samples = gsub("[a-zA-Z_.,]",'',colnames(counts))
	SNPgroups = gsub("[A-Za-z_.,[:punct:]]",'',geno[i,1])
	tgroup1 = which(samples %in% unlist(strsplit(SNPgroups,'\\s')))
	tgroup2 = which(samples %in% unlist(strsplit(SNPgroups,'\\s'))==FALSE)
	tgroups[tgroup1] = "ALT"
	tgroups[tgroup2] = "REF"
	tdesign=model.matrix(~tgroups)
	tcounts = DGEList(counts=tcounts,group=tgroups)
	tcounts$'samples'[,'lib.size'] = all_counts$'samples'[,'lib.size']
	tcounts$'samples'[,'norm.factors'] = all_counts$'samples'[,'norm.factors']
	tcounts = estimateGLMCommonDisp(tcounts,tdesign)
	tcounts = estimateGLMTrendedDisp(tcounts,tdesign)
	tcounts = estimateGLMTagwiseDisp(tcounts,tdesign)
	if (tcounts$'common.dispersion' > 0 & tcounts$'trended.dispersion' > 0 & tcounts$'tagwise.dispersion' > 0 ){
	tfit = glmFit(tcounts,tdesign)
	tlrt = glmLRT(tfit)
	tdisp = tlrt$'dispersion'
	# NULL MODEL
	rdesign=tdesign
	rdesign[,2]=1
	rcounts = DGEList(counts=tcounts,group=factor(rep(1,times=66)))
	rcounts$'samples'[,'lib.size'] = all_counts$'samples'[,'lib.size']
	rcounts = estimateGLMCommonDisp(rcounts,design=NULL)
	rcounts = estimateGLMTrendedDisp(rcounts,design=NULL)
	rcounts = estimateGLMTagwiseDisp(rcounts,design=NULL)
	rfit = glmFit(rcounts,design=NULL)	
	rfit$'design' = rdesign[,]
	rownames(rfit$'design') = colnames(rfit$'counts')	
	rlrt = glmLRT(rfit,coef=1)
	rdisp = rlrt$'dispersion'
	dispRatio = 1-sqrt(tdisp)/sqrt(rdisp)
	R2 = 1-exp(1)^((tlrt$'table'[i,3])/(-66))
	return(c(tlrt$'table',tdisp,rdisp,dispRatio,R2))
}
}

results=list()
for (i in genes){
	results[i] = tryCatch.W.E(
	lrt_test(i)$'value'
	)
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


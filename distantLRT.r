
'''
LRT test for distant effects. All expressed genes (with at least one count per 90% of samples) will be tested against all SNP genotypes.  
'''

snpPath = # PATH TO haplotype_parser_v3.py OUTPUT
haploPath = # PATH TO snp_parser_v2.py OUTPUT

'''
Lines (SNPs) from snp file and haplotype file will be parsed to select a non-redundant set of SNP genotypes. Preference will be give for genotypes as defined in haplotype file.
'''

outPath = # PATH TO DIRECTORY WHERE RESULTS WILL BE SAVED - ONE FILE WILL BE GENERATED PER FOCAL GENE

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
counts = DGEList(counts=counts)
counts = calcNormFactors(counts,method='TMM')

snp = read.table(snpPath,sep='\t',header=FALSE) 
haplotypes = read.table(haploPath,sep='\t')
rownames(haplotypes) = haplotypes[,1]

get_samples = function(x){
         return(unlist(strsplit(unlist(strsplit(as.character(x),'],'))[2],'\\['))[2])
        }

genotypes = unique(c(as.character(snp[,1]),as.character(haplotypes[,1])))
geno = matrix(nrow=length(genotypes),ncol=1)
rownames(geno) = genotypes
for (i in rownames(geno)){ # parsing of samples that represent two alternative maternal alleles
        if (i %in% rownames(haplotypes)){
                geno[i,1] = as.character(get_samples(haplotypes[i,7]))
        }
        else {
                if (i %in% snp[,1]){
                        geno[i,1] = as.character(snp[,6][snp[,1] %in% i][1])

                }
        }
}


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


lrt_test_trans = function(target,trans){
	tcounts = counts[target,]
	tgroups = {}
	samples = gsub("[a-zA-Z_.,]",'',colnames(counts))
	SNPgroups = gsub("[a-zA-Z_.,[:punct:]]",'',geno[trans,1])
	tgroup1 = which(samples %in% unlist(strsplit(SNPgroups,'\\s')))
	tgroup2 = which(samples %in% unlist(strsplit(SNPgroups,'\\s'))==FALSE)
	tgroups[tgroup1] = "G1"
	tgroups[tgroup2] = "G2"
	tdesign=model.matrix(~tgroups)
	tcounts = DGEList(counts=tcounts,group=tgroups)
	tcounts$'samples'[,'lib.size'] = counts$'samples'[,'lib.size']
	tcounts$'samples'[,'norm.factors'] = counts$'samples'[,'norm.factors']
	tcounts = estimateGLMCommonDisp(tcounts,tdesign)
	tcounts = estimateGLMTagwiseDisp(tcounts,tdesign)
	if (tcounts$'common.dispersion' > 0 & tcounts$'tagwise.dispersion' > 0 ){
		tfit = glmFit(tcounts,tdesign)
		tlrt = glmLRT(tfit)
		tdisp = tlrt$'dispersion'
		# NULL MODEL
		rdesign=tdesign
		rdesign[,2]=1
		rcounts = DGEList(counts=tcounts,group=factor(rep(1,times=66)))
		rcounts$'samples'[,'lib.size'] = counts$'samples'[,'lib.size']
		rcounts = estimateGLMCommonDisp(rcounts,design=NULL)
		rcounts = estimateGLMTagwiseDisp(rcounts,design=NULL)
		rfit = glmFit(rcounts,design=NULL)	
		rfit$'design' = rdesign[,]
		rownames(rfit$'design') = colnames(rfit$'counts')	
		rlrt = glmLRT(rfit,coef=1)
		rdisp = rlrt$'dispersion'
		dispRatio = 1-sqrt(tdisp)/sqrt(rdisp)
		R2 = 1-exp(1)^((tlrt$'table'[target,3])/(-66))
		return(list(tlrt$'table',tdisp,rdisp,dispRatio,R2))
	}
}


for (target in rownames(counts)){ # for each gene
	results=list()
	if (length(which(c(counts[target,] > 0))) > 0.9 * ncol(counts)){
		for (trans in rownames(geno)){ # for each trans gene
			if (length(which(c(counts$counts[trans,] > 0)))  > 0.9 * ncol(counts)){
				results[[target]][[trans]] = tryCatch.W.E(lrt_test_trans(target,trans))$'value'
			} 
		}
	}
	if (length(results[[target]]) > 0){
		resMat = matrix(nrow=length(results[[target]]),ncol = 8)
		resMat = as.data.frame(resMat)
		for (y in 1:length(results[[target]])){
			if (grepl(target,names(results[[target]])[y])==F & length(unlist(results[[target]][y])) == 8){
				resMat[y,1:8] = unlist(results[[target]][y])
				rownames(resMat)[y] = names(results[[target]])[y]
			}
		}		
		colnames(resMat)[1:8]=c('logFC','logCPM','LR','PValue','model_dispersion','null_dispersion','dispRatio','pseudo-R2')
		write.table(resMat,paste(outPath,target,sep=''))
	}
}


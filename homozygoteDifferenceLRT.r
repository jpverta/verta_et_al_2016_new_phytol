# R version 2.15.2

'''
Test of homozygote difference of haploid associations.
'''

genoPath = # PATH TO OUTPUT FROM genotypeFrequencies.r

countPath = # PATH TO COUNT DATA (output from htseq_count_summary.r), alternatively use the code below to summarize counts across HTSeq output files (.htseq_counts)

outPath = # PATH TO OUTFILE

library(stringr)
library(edgeR)
library(qvalue)

genoTable = read.table(genoPath,header=T)

test = ''' 

A VECTOR OF ASSOCIATIONS TO BE TESTED, FOCAL and EFFECT GENE NAMES SEPARATED WITH A SEMICOLON, i.e. focal_gene:effect_gene, 
IN THE CASE OF LOCAL TRANS EFFECTS THE TWO IDENTIFIERS ARE IDENTICAL, i.e. focal_gene:focal_gene

These effect genes must be found in .genot files, otherwise they will not be tested. 

'''

############### REVISED VERSION #################

haploList = read.table('haplotype_lists.out',sep='\t')
haploList = unique(haploList)
rownames(haploList) = haploList[,1]



HRsamples=list()
HAsamples=list()
HZsamples=list()
samples={}
for (i in list.files()){
    if (grepl('.genot',i)==T & grepl('NEW',i)==F & grepl('snp',i) == F & grepl('test',i) == FALSE & grepl('NM',i) == F & grepl('.txt',i) == F){
        print(i)
        geno = read.table(i,header=T,sep='\t')
        rownames(geno) = paste(geno[,1],geno[,2],sep=':')
        samples = append(samples,unlist(strsplit(i,'\\.'))[1])
        for (y in rownames(d)){
            gene = unlist(strsplit(y,':'))[1]
            effect = unlist(strsplit(y,':'))[2]
            if (length(which(genoTable[genoTable[,1] %in% effect,'Heteroz_obs']>=25))> 0.5*length(genoTable[genoTable[,1] %in% effect,'Heteroz_obs']) &
            length(which(genoTable[genoTable[,1] %in% effect,'Heteroz_obs']<=41))>0.5* length(genoTable[genoTable[,1] %in% effect,'Heteroz_obs'])){
                if (length(which(geno[grepl(effect,rownames(geno)),'Embryo_AF'] == '[0.5]')) > length(geno[grepl(effect,rownames(geno)),'Embryo_AF'] == '[0.5]')/2){
                    HZsamples[[y]] = append(HZsamples[[y]],unlist(strsplit(i,'\\.'))[1])
                }
                if (length(which(geno[grepl(effect,rownames(geno)),'Embryo_AF'] == '[1.0]')) > length(geno[grepl(effect,rownames(geno)),'Embryo_AF'] == '[1.0]')/2){
                    g = geno[grepl(effect,rownames(geno)),]
                    
                    haploA = unlist(strsplit(unlist(strsplit(str_replace_all(haploList[effect,2],"[[:punct:]]",""),',')),'\\s'))
                    haploB = unlist(strsplit(unlist(strsplit(str_replace_all(haploList[effect,3],"[[:punct:]]",""),',')),'\\s'))
                    
                    # GENOTYPE CALLING FUNCTION
                    if (length(haploB) > 0){
                        # require that 90% of the SNP calls within a set of linked SNPs are belongin to the same haplotype
                        if (length(which(g[g[,2] %in% haploA,4] == g[g[,2] %in% haploA,5])) > 0.9* length(haploA) &
                        length(which(g[g[,2] %in% haploB,3] == g[g[,2] %in% haploB,5])) > 0.9*length(haploB)){
                            HAsamples[[y]] = append(HAsamples[[y]],unlist(strsplit(i,'\\.'))[1])
                        }
                        else {
                            if (length(which(g[g[,2] %in% haploA,3] == g[g[,2] %in% haploA,5])) > 0.9* length(haploA) &
                            length(which(g[g[,2] %in% haploB,4] == g[g[,2] %in% haploB,5])) > 0.9*length(haploB)){
                                HRsamples[[y]] = append(HRsamples[[y]],unlist(strsplit(i,'\\.'))[1])
                            }
                        }
                    }
                    else{
                        if (length(which(g[g[,2] %in% haploA,4] == g[g[,2] %in% haploA,5])) > 0.9*length(haploA)){
                            HAsamples[[y]] = append(HAsamples[[y]],unlist(strsplit(i,'\\.'))[1])
                        }
                        else {
                            if (length(which(g[g[,2] %in% haploA,3] == g[g[,2] %in% haploA,5])) > 0.9*length(haploA)){
                                HRsamples[[y]] = append(HRsamples[[y]],unlist(strsplit(i,'\\.'))[1])
                            }
                        }
                        
                    }
                }
            }
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


counts = read.table(countPath)

cnames = {}
for (i in colnames(counts)){
        if (grepl('Index',i)==FALSE){
                cnames=append(cnames,str_replace_all(i,"[a-z_.,[:punct:]]",""))
        }
        else{
                cnames = append(cnames,unlist(strsplit(i,'\\.'))[5])
}
}

cnames = str_replace_all(cnames,"[A-Z]","")
colnames(counts) = cnames
counts = counts[samples]

all_counts = DGEList(counts=counts)
all_counts = calcNormFactors(all_counts,method='TMM')

lrt_trans_test = function(i){
        samples = gsub("[a-zA-Z_.,]",'',colnames(counts))
        gene = unlist(strsplit(i,':'))[1]
        tcounts1 = counts[gene,which(samples %in% unlist(HAsamples[as.character(i)]))]
        tcounts2 = counts[gene,which(samples %in% unlist(HRsamples[as.character(i)]))]
        tcounts = c(as.numeric(tcounts1),as.numeric(tcounts2))
        tcounts = t(as.matrix(tcounts))
        tgroups = {}
        tgroup1 = rep('A',times=length(tcounts1))
        tgroup2 = rep('B',times=length(tcounts2))
        tgroups = c(tgroup1,tgroup2)
        tdesign=model.matrix(~tgroups)
        tcounts = DGEList(counts=tcounts,group=tgroups)
		tcounts$'samples'[,'lib.size'] = all_counts$'samples'[c(which(rownames(all_counts$'samples') %in% unlist(HAsamples[as.character(i)])),which(rownames(all_counts$'samples') %in% unlist(HRsamples[as.character(i)]))),'lib.size']
		tcounts$'samples'[,'norm.factors'] = all_counts$'samples'[c(which(rownames(all_counts$'samples') %in% unlist(HAsamples[as.character(i)])),which(rownames(all_counts$'samples') %in% unlist(HRsamples[as.character(i)]))),'norm.factors']
        tcounts = estimateGLMCommonDisp(tcounts,tdesign)
        tcounts = estimateGLMTagwiseDisp(tcounts,tdesign)
        tfit = glmFit(tcounts,tdesign)
        et = glmLRT(tfit)
return(c(mean(as.numeric(tcounts1)),mean(as.numeric(tcounts2)),length(unlist(HAsamples[i])),length(unlist(HRsamples[i])),length(unlist(HZsamples[i])),et$'table'$'logFC',et$'table'$"PValue"))
}


res=c()
homozTable = data.frame()
for (i in test){
        if (length(unlist(HRsamples[as.character(i)])) <= 25 & length(unlist(HRsamples[as.character(i)])) >= 10
        & length(unlist(HAsamples[as.character(i)])) <= 25 & length(unlist(HAsamples[as.character(i)])) >= 10){
        res = tryCatch.W.E(lrt_trans_test(i))$'value'
        if (length(res) == 7){
        homozTable[i,1] = res[1]
        homozTable[i,2] = res[2]
        homozTable[i,3] = res[3]
        homozTable[i,4] = res[4]
		homozTable[i,5] = res[5]
        homozTable[i,6] = res[6]
		homozTable[i,7] = res[7]
}
}
}

write.table(homozTable,outPath)



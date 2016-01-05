
'''
Output parsing and plots for local (1n) - cis (2n) comparison.
'''

localPath =  '/fml/jones/home/jverta/PhD/analysis/megagametophyte/local_LRT_htseq_counts.txt' # output from localLRT.r

l2nPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffLocalHaploCounts.txt' # output from hmoozygoteDifferenceLRT.r for local effects

#tissueCompPath =  # Output from tissueComparison.r

asePath = '/fml/jones/home/jverta/PhD/analysis/embryo/aseGLMtable_corrected_all_closeDrop.txt' # Output from haploASE.r or singlesASE.r (or a combined table of the two outputs)
# aseGLMtable_corrected_all_closeDrop.txt

library(qvalue)

#REVISED - INCLUDE NON-PREFERENTIALLY EXPRESSED GENES!
#non-preferentially expressed genes
#tissue = read.table(tissueCompPath)
#adj = qvalue(tissue[,4])
#sigT = tissue[adj$qvalue<0.01,]
#pref = sigT[sigT[,1] < -2 & sigT[,1] < 0 | sigT[,1] > 0 & sigT[,1] > 2,]
#np = rownames(tissue)[rownames(tissue) %in% rownames(pref)==F]

# test tables

ase = read.table(asePath)
rownames(ase) = ase[,1]
adjA = qvalue(ase[,7])
sigA = ase[adjA$qvalue<0.01,]

local = read.table(localPath)
local = local[is.na(local[,4])==F,]
adjL = qvalue(local[,'PValue'])
sigL = local[adjL$qvalue<0.01,]

l2N = read.table(l2nPath)
l2N = l2N[l2N[,3] > 9 & l2N[,3] < 25 & l2N[,4] > 9 & l2N[,4] < 25,]
effects = c()
genes=c()
for (i in rownames(l2N)){
	g = unlist(strsplit(i,':'))[1]
	e = unlist(strsplit(i,':'))[2]
	effects = append(effects,e)
	genes = append(genes,g)
}
rownames(l2N) = genes
adjL2N = qvalue(l2N[,'PValue'])
sigL2N = l2N[adjL2N$qvalue<0.01,]

# "local" == local effects tested in 1N
# "sigL" == significant local 1N effects
# "ase" == all tested cis effects
# "sigA" == significant cis effects
# "local2N" == local effects tested in 2N
# "sigL2N" == significant local 2N effects

## LOCAL TRANS: difference between homozygous 2n but no differences in heterozygous 2n
# filter for genes that change signs between 1N local and 2N homozygote differences

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# genes tested for all three effects, local, local 2N and cis
ml = rownames(local)[rownames(local) %in% rownames(l2N) & rownames(local) %in% ase[,1]]

# sign in genes under local effects in 1N and 2N homozygotes -> these need to be the same

l1s = local[ml,'logFC'] > 0
l2s = l2N[ml,'logFC'] > 0

length(which(l1s & l2s | l1s ==F & l2s ==F))
length(ml)

samesignL = ml[which(l1s & l2s | l1s ==F & l2s ==F)]

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# genes under local trans effects: local effects in 1N and 2N but no cis effects and tested for all
length(rownames(sigL)[rownames(sigL) %in% sigA[,1]==F & rownames(sigL) %in% ase[,1] & rownames(sigL) %in% rownames(sigL2N) ]) #& rownames(sigL) %in% samesignL])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# tested genes:
length(rownames(sigL)[rownames(sigL) %in% ase[,1] & rownames(sigL) %in% rownames(l2N) ]) #& rownames(sigL) %in% samesignL])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# genes under local trans effects and local effects observed in all three genotypes
length(rownames(sigL2N)[rownames(sigL2N) %in% sigA[,1]==F ]) #& rownames(sigL2N) %in% samesignL])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# tested genes
length(rownames(sigL2N)[rownames(sigL2N) %in% ase[,1] ]) #& rownames(sigL2N) %in% samesignL])

# LOCAL CIS: difference between local 1n and difference in heterozygous 2n

# genes tested for local and cis
mc = rownames(local)[rownames(local) %in% ase[,1]]

# sign in genes under local effects in 1N and 2N heterozygotes -> these need to be the same

l1s = local[mc,1] > 0
l2s = ase[mc,6] > 0

samesignA = mc[which(l1s & l2s | l1s ==F & l2s ==F)]

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# local cis
length(rownames(sigL)[rownames(sigL) %in% rownames(sigA)]) #& rownames(sigL) %in% samesignA])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# local tested for cis
length(rownames(sigL)[rownames(sigL) %in% rownames(ase) ]) #& rownames(sigL) %in% samesignA])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# cis local
length(rownames(sigA)[rownames(sigA) %in% rownames(sigL) ]) #& rownames(sigA) %in% samesignA])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# cis tested for local 
length(rownames(sigA)[rownames(sigA) %in% rownames(local) ]) #& rownames(sigA) %in% samesignA])

# LOCAL HOMOZ. DIFFERENCE

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# genes tested for local, local 2N
mh = rownames(local)[rownames(local) %in% rownames(l2N)]

# sign in genes under local effects in 1N and 2N homozygotes -> these need to be the same

l1s = local[mh,1] > 0
l2s = l2N[mh,6] > 0

samesignH = mh[which(l1s & l2s | l1s ==F & l2s ==F)]

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# local homoz. difference
length(rownames(sigL)[rownames(sigL) %in% rownames(sigL2N) ]) #& rownames(sigL) %in% samesignH])

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# local tested for homoz. difference
length(rownames(sigL)[rownames(sigL) %in% rownames(l2N) ]) #& rownames(sigL) %in% samesignH])





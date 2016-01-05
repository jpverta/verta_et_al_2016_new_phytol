

'''
Output parsing and plots for homozygote difference (2n) - cis/trans (2n) comparison.
'''

#tissueCompPath =  # Output from tissueComparison.r

lgPath = '/fml/jones/home/jverta/PhD/analysis/megagametophyte/LG_RF035_LOD6.txt' # Output from mapFiles.r (linkage blocks)

#localPath1 =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffSNPcountsCloseDrop_1half_local.txt' # Output from homozygoteDifferenceLRT.r - SNP counts
#localPath2 =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffSNPcountsCloseDrop_2half_local.txt' # Output from homozygoteDifferenceLRT.r - SNP counts

localPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffLocalHaploCounts.txt' # Output from homozygoteDifferenceLRT.r - Haplotype counts

#localPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffSNPcountsCloseDropLocal_TEST.txt' # Output from homozygoteDifferenceLRT.r - SNPs/haplotype counts

#localPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffSNPcountsLocal.txt'
 
#localPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffSNPcountsLocal_SNPonly.txt'

# REDO THIS FOR HOMOZYGOTE DIFFERENCE CALCULATED BASED ON SNP COUNTS IN DIPLOID HOMOZYGOTES!

asePath =  '/fml/jones/home/jverta/PhD/analysis/embryo/aseGLMtable_corrected_all_closeDrop.txt' # Output from haploASE.r or singlesASE.r (or a combined table of the two outputs)

# ONLY haplotype-based ASE counts
#asePath =  '/fml/jones/home/jverta/PhD/analysis/embryo/haploASEglm_closeDrop.txt' # Output from haploASEgml_closeDrop.r

haploPath = '/fml/jones/home/jverta/PhD/analysis/embryo/' # PATH TO HAPLOTYPE TABLE (output from haplotypeHeterozygoteFrequencies.r)

library(qvalue)

# REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL GENES
#non-preferentially expressed genes
#tissue = read.table(tissueCompPath)
#adj = qvalue(tissue[,4])
#sigT = tissue[adj$qvalue<0.01,]
#pref = sigT[sigT[,1] < -2 & sigT[,1] < 0 | sigT[,1] > 0 & sigT[,1] > 2,]
#np = rownames(tissue)[rownames(tissue) %in% rownames(pref)==F]
#length(np)

# test tables for homozygote difference and cis effects
ase = read.table(asePath,stringsAsFactors=F)
rownames(ase) = ase[,1]
#local1 = read.table(localPath1)
#local2 = read.table(localPath2)
#colnames(local1) = colnames(local2)
#local = rbind(local1,local2)
local = read.table(localPath,stringsAsFactors=F)
local = local[local[,3] > 9 & local[,3] < 25 & local[,4] > 9 & local[,4] < 25,]

# If considering table for gene counts:
#colnames(local) = c('AO','BO','HAsamples','HBsamples','HZsamples','logFC','Pvalue')

# REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL GENES
#effects = c()
#genes=c()
#for (i in rownames(local)){
#	g = unlist(strsplit(i,':'))[1]
#	e = unlist(strsplit(i,':'))[2]
#	effects = append(effects,e)
#	genes = append(genes,g)
#}
#rownames(local) = genes
#local=local[rownames(local) %in% np,]

adjL = qvalue(local[,'PValue']) # SNP counts
sigL = local[adjL$qvalue<0.01,]

# HAPLO ASE
#ase = ase[is.na(ase[,'LRT_pvalue'])==F,]
#adjA = qvalue(ase[,'LRT_pvalue'])

ase = ase[is.na(ase[,'exactTest_pvalue'])==F,]
adjA = qvalue(ase[,'exactTest_pvalue'])
sigA = ase[adjA$qvalue<0.01,]

# "local" == local effects tested in 2N
# "sigL" == significant local 2N effects
# "ase" == all tested cis effects
# "sigA" == significant cis effects

## signs in 2N heterozygotes and homozygotes

# REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL GENES
m = rownames(ase)[rownames(ase) %in% rownames(local)]
#m = as.character(ase[ase[,1] %in% rownames(local),1])
#m = m[m %in% np]

signA = ase[m,'logFC']<0
signL = local[m,'logFC']<0

signMatch = m[signA ==T & signL ==T |  signA ==F & signL ==F]
signSwitch = m[m %in% signMatch ==F]

signMatch
signSwitch

write.table(signMatch,'localCisSameSignHaplo2NCounts.txt')

# overlap

length(rownames(sigL)[rownames(sigL) %in% sigA[,1] & rownames(sigL) %in% signMatch]) # genes that show homozygote difference and cis effects

length(rownames(sigL)[rownames(sigL) %in% sigA[,1] == F & rownames(sigL) %in% ase[,1]])



############### COMPENSATORY ###################

# compensatory -> cis effects but no homozygote difference, tested for homozygote difference

# REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL GENES
tg = as.character(sigA[sigA[,1] %in% rownames(local) & sigA[,1] %in% rownames(sigL) == F,1])# & sigA[,1] %in% signMatch,1])
length(tg) # compensatory

# REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL GENES
tested = as.character(sigA[sigA[,1] %in% rownames(local),1])# & sigA[,1] %in% signMatch,1])
length(tested) # number of tested genes

pdf('compensatory_haploCounts_annotated.pdf')
plot(local[tg,'logFC'], sigA[tg,'logFC'],xlim=c(-2,2),ylim=c(-2,2),pch=16,cex.lab=2,cex.axis=2)
abline(h=0,v=0)
abline(0,1,col='darkgrey',lty=2,lwd=2.5)
#dev.off()

# mark points with P<0.05 without qvalue correction in homozygote comparison
points(local[tg[local[tg,'PValue'] < 0.05],'logFC'], sigA[tg[local[tg,'PValue'] < 0.05],'logFC'],col='green',pch=16)

#points(local[tg[local[tg,'logFC'] > sigA[tg,'logFC'] & local[tg,'logFC'] > 0 & sigA[tg,'logFC'] > 0],'logFC'],sigA[tg[local[tg,'logFC'] > sigA[tg,'logFC'] & local[tg,'logFC'] > 0 & sigA[tg,'logFC'] > 0],'logFC'],col='red',pch=16)

#points(local[tg[local[tg,'logFC'] < sigA[tg,'logFC'] & local[tg,'logFC'] < 0 & sigA[tg,'logFC'] < 0],'logFC'],sigA[tg[local[tg,'logFC'] < sigA[tg,'logFC'] & local[tg,'logFC'] < 0 & sigA[tg,'logFC'] < 0],'logFC'],col='blue',pch=16)

dev.off()

# pseudo R2 in contradicting log FC cases
ase[tg[local[tg,'logFC'] > sigA[tg,'logFC'] & local[tg,'logFC'] > 0 & sigA[tg,'logFC'] > 0],'pseudo.R2']
ase[tg[local[tg,'logFC'] < sigA[tg,'logFC'] & local[tg,'logFC'] < 0 & sigA[tg,'logFC'] < 0],'pseudo.R2']



# CIS
cis = as.character(sigA[sigA[,1] %in% rownames(local) & sigA[,1] %in% rownames(sigL),1])
x11()
plot(local[cis,'logFC'], sigA[cis,'logFC'],xlim=c(-2,2),ylim=c(-2,2),pch=16,cex.lab=2,cex.axis=2)
abline(h=0,v=0)
abline(0,1,col='darkgrey',lty=2,lwd=2.5)

### plot of all three categories

# "local" == local effects tested in 2N
# "ase" == all tested cis effects
# "sigA" == significant cis effects

#compensatory (cis effect (2N) but no local effect (1N))
rownames(sigA)=sigA[,1]
tg = as.character(sigA[sigA[,1] %in% rownames(local) & sigA[,1] %in% rownames(sigL) == F,1])
#local trans
lt = rownames(sigL[rownames(sigL) %in% rownames(ase) & rownames(sigL) %in% sigA[,1]==F,])
#local cis
lc = rownames(sigL[rownames(sigL) %in% sigA[,1],])

##
write.table(lc,'localCisHaploCounts.txt')
write.table(lt,'localTransHaploCounts.txt')
write.table(tg,'localCompHaploCounts.txt')
##


library(ggplot2)

FCplot = data.frame(local[,'logFC'],ase[rownames(local),'logFC'])
rownames(FCplot) = rownames(local)

FCplot[rownames(FCplot) %in% tg,'class'] = 'tg'
FCplot[rownames(FCplot) %in% tg,'shape'] = 1
FCplot[rownames(FCplot) %in% lc,'class'] = 'lc'
FCplot[rownames(FCplot) %in% lc,'shape'] = 0
FCplot[rownames(FCplot) %in% lt,'class'] = 'lt'
FCplot[rownames(FCplot) %in% lt,'shape'] = 2
FCplot[is.na(FCplot[,'class']),'class'] = 'none'
FCplot[FCplot[,'class'] == 'none','shape'] = 4

ggplotColours <- function(n=6, h=c(0, 360) +15){
    if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
    hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

ggplotColours(n=3)
# "#F8766D" "#00BA38" "#619CFF"

#ggplot(FCplot,aes(local....logFC..,ase.rownames.local....logFC..)) + geom_point(aes(colour = class,fill=NA,shape=factor(shape)))

ggplot(FCplot,aes(local....logFC..,ase.rownames.local....logFC..)) + geom_point(aes(colour=class,shape = factor(shape)), fill='NA') + scale_color_manual(labels = c('lc','lt','none','tg'),values = c('#F8766D','#00BA38','grey','#619CFF')) #+ scale_shape_manual(values = FCplot$shape)


#pdf('local_vs_cis_haploCounts.pdf')
#pdf('local_vs_cis_SNPcounts.pdf')
ggplot(FCplot,aes(local....logFC..,ase.rownames.local....logFC..)) + geom_point(data = FCplot[FCplot$class == 'none',],colour = '#999999',shape=4) + geom_point(data = FCplot[FCplot$class %in% 'none' ==F,], aes(colour=class,shape = factor(shape)), fill='NA') + scale_color_manual(values = c('#F8766D','#00BA38','#619CFF')) + theme(legend.position = "none") + ylim(-3,3) + xlim(-3,3)
#dev.off()

	
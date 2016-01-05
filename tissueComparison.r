
'''
LRT comparison of megagametophyte and embryo expression. 
'''

outPath = # PATH TO OUTFILE

mcountPath = # PATH TO MEGAGAMETOPHYTE COUNT DATA (output from htseq_count_summary.r), alternatively use the code below to summarize counts across HTSeq output files (.htseq_counts)

ecountPath = # PATH TO EMBRYO COUNT DATA (output from htseq_count_summary.r), alternatively use the code below to summarize counts across HTSeq output files (.htseq_counts)



haploid = read.table('mgg_htseq_counts.txt')
diploid = read.table('embryo_htseq_counts.txt')
diploid = diploid[,1:66]


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

mcounts = read.table(mcountPath,header=T)
ecounts = read.table(ecountPath,header=T)
tcounts = cbind(mcounts,ecounts)
tcounts = DGEList(counts=tcounts)
tcounts = calcNormFactors(tcounts,method='TMM')

tgroups = c(rep("M",times=ncol(mcounts)),rep('E',times=ncol(ecounts)))
tdesign = model.matrix(~tgroups)
tcounts = estimateCommonDisp(tcounts,tdesign)
tcounts = estimateTagwiseDisp(tcounts)
tfit = glmFit(tcounts,tdesign)
tlrt = glmLRT(tfit)

write.table(tlrt$table,outPath)

# difference between gene groups --> not working at the moment 

compgenes = read.table('localCompHaploCounts.txt')
cisgenes = read.table('localCisHaploCounts.txt')
transgenes = read.table('localTransHaploCounts.txt')

#subset haploid table to categories
#rbind categories together
#draw up design table
#estimate dispersions
#test DE

hapComp = haploid[rownames(haploid) %in% compgenes[,1],]
hapCis = haploid[rownames(haploid) %in% cisgenes[,1],]
hapTrans = haploid[rownames(haploid) %in% transgenes[,1],]

tcounts = rbind(hapComp,hapCis,hapTrans)
tgroups = c(rep("Comp",times=nrow(hapComp)),rep("Cis",times=nrow(hapCis)),rep("Trans",times = nrow(hapTrans)))
tcounts = DGEList(counts=t(tcounts))
tdesign = model.matrix(~tgroups)
tcounts = estimateCommonDisp(tcounts,tdesign)
tcounts = estimateTagwiseDisp(tcounts)
tfit = glmFit(tcounts,tdesign)
tlrt = glmLRT(tfit)

res = tlrt$table

qval = qvalue(res[,'PValue'])


# RPKM plots
library(ggplot2)

diploid = read.table('embryo_rpkm_counts.txt')
haploid = read.table('mgg_rpkm_counts.txt')

ge = data.frame(apply(haploid,1,mean), apply(diploid,1,mean))
ge = ge[ge[,1] > 0 & ge[,2] > 0,]

pdf('haploid_vs_diploid_rpkm.pdf')
ggplot(ge, aes(ge[,1],ge[,2]))+geom_point()+geom_smooth(method=lm)+theme(text = element_text(size=20))+xlab('Haploid RPKM')+ylab('Diploid RPKM')
dev.off()




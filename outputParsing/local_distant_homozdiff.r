'''
Output parsing and plot for local and distant effects (1n) vs. homozygote difference (2n)
'''

distPath =  '/fml/jones/home/jverta/PhD/analysis/megagametophyte/trans_LRT_significative_rf_combBlocks_R2.txt' # Output from distantLRTsummary.r, significan effects

d2nSigPath =  # output from homozygoteDifferenceLRT.r for distant effects, significant effects

d2nPath =  '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffDistantHTseqCounts_withSNP_combBlocks.txt' # output from homozygoteDifferenceLRT.r for distant effects

localPath =  '/fml/jones/home/jverta/PhD/analysis/megagametophyte/local_LRT_htseq_counts.txt' # output from localLRT.r

l2nPath = '/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffLocalGenesHTseqCounts.txt' # output from homozygoteDifferenceLRT.r for local effects

l1nPath =  '/fml/jones/home/jverta/PhD/analysis/megagametophyte/local_LRT_htseq_counts.txt' # Output from local_LRT.r

lgPath = '/fml/jones/home/jverta/PhD/analysis/megagametophyte/LG_RF035_LOD6.txt' # Output from mapFiles.r (linkage blocks)

#tissueCompPath =  # Output from tissueComparison.r

############################################

library(qvalue)

#REVISED - INCLUDE NON-PREFERENTIALLY EXPRESSED GENES!
#non-preferentially expressed genes
#tissue = read.table(tissueCompPath)
#adj = qvalue(tissue[,4])
#sigT = tissue[adj$qvalue<0.01,]
#pref = sigT[sigT[,1] < -2 & sigT[,1] < 0 | sigT[,1] > 0 & sigT[,1] > 2,]
#np = rownames(tissue)[rownames(tissue) %in% rownames(pref)==F]

lg = read.table(lgPath)

##### DISTANT VS. HOMOZYGOTE DIFFERENCE #####

dist = read.table(distPath,header=T,stringsAsFactors=F)
dist=unique(dist)
rownames(dist) = paste(dist[,1],dist[,2],sep=':')
dist[,11] = lg[dist[,1],2]
dist[,12] = lg[dist[,2],2]

alltrans=read.table(d2nPath,fill=T,row.names=NULL,stringsAsFactors=F)
alltrans = alltrans[as.numeric(alltrans[,4]) >9 & as.numeric(alltrans[,4]) < 25 & as.numeric(alltrans[,5]) >9 & as.numeric(alltrans[,5]) < 25,]
alltrans=unique(alltrans)
alltrans = alltrans[is.na(alltrans[,1])==F,]
rownames(alltrans) = alltrans[,1]
alltrans = alltrans[,-1]

effects = c()
genes=c()
for (i in rownames(alltrans)){
    g = unlist(strsplit(i,':'))[1]
    e = unlist(strsplit(i,':'))[2]
    effects = append(effects,e)
    genes = append(genes,g)
}
alltrans[,8] = lg[genes,2]
alltrans[,9] = lg[effects,2]


# significant 2N homoz effects

q = qvalue(alltrans[,7])
trans = alltrans[q$qvalue<0.01,]
nrow(trans)

# ALL TESTED EFFECTS

#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# overlap in gnee set tested for distant effects and 2N homoz. difference
mt = rownames(alltrans)[rownames(alltrans) %in% rownames(dist)]

ds = sign(dist[mt,'logFC'])
ts = sign(as.numeric(alltrans[mt,'logFC']))

samesign2N = mt[ds == ts]
length(samesign2N)

# tested effects influencing mapped focal loci and acting in same direction in 1N and 2N 
length(samesign2N[is.na(alltrans[samesign2N,8])])

# tested effects influencing non-mapped focal loci and acting in same direction in 1N and 2N 
length(samesign2N[is.na(alltrans[samesign2N,8])==F])

# ALL TESTED NON-REDUNDANT EFFECTS

effects = c()
genes=c()
for (i in rownames(alltrans)){
    g = unlist(strsplit(i,':'))[1]
    e = unlist(strsplit(i,':'))[2]
    effects = append(effects,e)
    genes = append(genes,g)
}
R=data.frame()
for (i in unique(genes)){
	dr = alltrans[which(genes %in% i),]
	if (nrow(dr) > length(unique(dr[,9]))){
		t=matrix(ncol=9,nrow=1)
		colnames(t) = colnames(dr)
		t[1,1:9] = NA
		dr = dr[order(as.numeric(dr[,7])),] # start with the lowest P-value associations
		for (y in 1:nrow(dr)){
			if (dr[y,9] %in% t[,9] == F){
				t = rbind(t,dr[y,])
			}
		}
		R = rbind(R,t)
	}
	else {
		R = rbind(R,dr)
	}
}	
R = R[is.na(R[,1])==F,]
nrow(R)

distHaploidR2 = dist[rownames(R),'pseudo.R2']


#REVISED - NOT FILTERING FOR TISSUE PREFERENTIAL EXPRESSION
# overlap in reduced set of  tested 2N homoz difference and distant effects
mr = rownames(R)[rownames(R) %in% rownames(dist)]

ds = sign(dist[mr,'logFC'])
ts = sign(as.numeric(R[mr,'logFC']))

samesign2NR = mr[ds == ts]

# tested focal non-mapped
length(samesign2NR[is.na(R[samesign2NR,8])])

# tested focal mapped
length(samesign2NR[is.na(R[samesign2NR,8])==F])


# SIGNIFICANT EFFECTS IN 2N

effects = c()
genes=c()
for (i in rownames(trans)){
    g = unlist(strsplit(i,':'))[1]
    e = unlist(strsplit(i,':'))[2]
    effects = append(effects,e)
    genes = append(genes,g)
}
R=data.frame()
for (i in unique(genes)){
	dr = trans[which(genes %in% i),]
	if (nrow(dr) > length(unique(dr[,9]))){
		t=matrix(ncol=9,nrow=1)
		colnames(t) = colnames(dr)
		t[1,1:9] = NA
		dr = dr[order(as.numeric(dr[,7])),] # start with the lowest P-value associations
		for (y in 1:nrow(dr)){
			if (dr[y,9] %in% t[,9] == F){
				t = rbind(t,dr[y,])
			}
		}
		R = rbind(R,t)
	}
	else {
		R = rbind(R,dr)
	}
}	
R = R[is.na(R[,1])==F,]
nrow(R)

distBothPloidiesR2 = dist[rownames(R),'pseudo.R2']

# R2 LR

t.test(distHaploidR2[distHaploidR2 %in% distBothPloidiesR2 ==F],distBothPloidiesR2)



#### LOCAL VS. HOMOZYGOTE DIFFERENCE ####
# REVISED FOR SNP-BASED 2N COUNTS
# local effects tested for 2N homoz. diff.
local = read.table(l2nPath)
local = local[local[,3] > 9 & local[,3] < 25 & local[,4] > 9 & local[,4] < 25,]

########################
## LOCAL HAPLO COUNTS ##
########################
local = read.table('/fml/jones/home/jverta/PhD/analysis/embryo/homozDiffLocalHaploCounts.txt')
local = local[local[,3] > 9 & local[,3] < 25 & local[,4] > 9 & local[,4] < 25,]
########################


#effects = c()
genes=c()
for (i in rownames(local)){
        g = unlist(strsplit(i,':'))[1]
        e = unlist(strsplit(i,':'))[2]
        effects = append(effects,e)
        genes = append(genes,g)
}
rownames(local) = genes
#local=local[rownames(local) %in% np,]

adjL = qvalue(local[,'PValue'])
sig2N = local[adjL$qvalue<0.01,]

# local effects
local1N = read.table(l1nPath)
local1N = local1N[is.na(local1N[,4])==F,]
adj= qvalue(local1N[,4])
sig1N = local1N[adj$qvalue<0.01,]

# local effects tested for 2N homoz. diff. 
lt = rownames(sig1N)[rownames(sig1N) %in% rownames(local)]

# sign

l1s = sign(local1N[lt,'logFC'])
l2s = sign(local[lt,'logFC'])

samesignL = lt[which(l1s == l2s)]

# number of tests
length(lt[lt %in% samesignL])

# 1N local effects that correspond to 2N homozygote difference
ml = rownames(sig1N)[rownames(sig1N) %in% rownames(sig2N)]

# significant local and 2N homoz.diff
length(ml[ml %in% samesignL])

# R2 LR

l1nPath =  '/fml/jones/home/jverta/PhD/analysis/megagametophyte/local_LRT_htseq_counts.txt' # Output from local_LRT.r
local = read.table(l1nPath)

localBothPloidiesR2 = local[ml,'pseudo.R2']

localHaploidR2 = local[lt[lt%in%ml==F],'pseudo.R2']

t.test(localHaploidR2,localBothPloidiesR2)

#write.table(samesignL,'local1N2NsameSignHaploCounts.txt')


#### PLOTS 

library(ggplot2)

d1NnomapFC = dist[rownames(dist) %in% mnomap & rownames(dist) %in% samesign,'logFC']
d1NmapFC = dist[rownames(dist) %in% mmap & rownames(dist) %in% samesign,'logFC']

d2NnomapFC = trans[rownames(trans) %in% mnomap & rownames(trans) %in% samesign,'logFC']
d2NmapFC = trans[rownames(trans) %in% mmap & rownames(trans) %in% samesign,'logFC']

l1NFC = sig1N[rownames(sig1N) %in% ml & rownames(sig1N) %in% samesignL,'logFC']
l2NFC = sig2N[rownames(sig2N) %in% ml & rownames(sig2N) %in% samesignL,'logFC']

plotFC = data.frame(haploid=c(d1NnomapFC,d1NmapFC,l1NFC),diploid=c(d2NnomapFC,d2NmapFC,l2NFC))

plotFC$class = c(rep('dnomap',times=length(d1NnomapFC)),rep('dmap',times=length(d1NmapFC)),rep('local',times=length(l1NFC)))

ggplot(plotFC,aes(haploid,diploid)) + geom_point(aes(colour=class))

pdf('local_distant_v_homozdiff_sameSign_closeUp.pdf')
ggplot(plotFC,aes(haploid,diploid)) + geom_point(data = plotFC[plotFC$class == 'dnomap',],aes(colour = '#F8766D',fill='NA'),shape=1)  + geom_point(data = plotFC[plotFC$class == 'local',], aes(colour='#00BA38',fill='NA'),shape=2) + geom_point(data = plotFC[plotFC$class == 'dmap',], aes(colour='#619CFF',fill='NA'),shape=0) + theme(text = element_text(size=22)) + theme(legend.position = "none")
dev.off()

pdf('local_distant_v_homozdiff_sameSign.pdf')
ggplot(plotFC,aes(haploid,diploid)) + geom_point(data = plotFC[plotFC$class == 'dnomap',],aes(colour = '#F8766D',fill='NA'),shape=1)  + geom_point(data = plotFC[plotFC$class == 'local',], aes(colour='#00BA38',fill='NA'),shape=2) + geom_point(data = plotFC[plotFC$class == 'dmap',], aes(colour='#619CFF',fill='NA'),shape=0) + theme(text = element_text(size=22)) + ylim(-2,2) + xlim(-2,2) + theme(legend.position = "none")
dev.off()

# local

x11()
plot(sig1N[rownames(sig1N) %in% rownames(sig2N),'logFC'],sig2N[rownames(sig2N) %in% rownames(sig1N),'logFC'],pch=23,col='red',cex=0.7)
points(sig1N[rownames(sig1N) %in% ml & rownames(sig1N) %in% samesignL,1] ==F ,sig2N[rownames(sig2N) %in% ml & rownames(sig2N) %in% samesignL,6] ==F ,pch=23,col='grey',cex=0.7)


#pdf('local_distant_v_homozdiff_sameSign.pdf')

#dev.off()





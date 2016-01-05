'''
Matrix plot of local and distant effects.
'''

library(qvalue)

localPath = '/Volumes/LaCie/Data/Katak/Megagametophyte/eQTL/local_LRT_htseq_counts.txt' # output from localLRT.r
distPath = '/Volumes/LaCie/Data/Katak/Megagametophyte/eQTL/trans_LRT_significative_rf_combBlocks_R2.txt' # output from distantLRTsummary.r
mapPath = '/Volumes/LaCie/Data/Katak/Megagametophyte/eQTL/LG_RF035_LOD6.txt' # output from mapFiles.r (linkage blocks)
rfPath = '/Users/jverta/Documents/PhD/analysis/recombination_fractions_rf035_lod6.txt' # output from mapFiles.r (recombination fractions)


local = read.table(localPath)
local = local[is.na(local[,4]) == F,]
adj = qvalue(local[,4])
loc=local[adj$qvalue<0.01,]

dist = read.table(distPath,header=T,stringsAsFactors=F)

lg = read.table(mapPath)

rf = read.table(rfPath)

# select only mapped pairs 

dist_map = dist[is.na(dist[,7]) ==F,]
genes = as.character(dist_map[,1])
effects = as.character(dist_map[,2])

#############
## ggplot2 ##
#############

# order genes within lg's according to recombination fraction

lgOrdered=data.frame()
for (i in c(1:max(lg[,2]))){
    genes = lg[lg[,2]==i,]
    r = rf[rownames(genes),]
    genes = genes[order(r[,1],r[,2],r[,3]),]
    lgOrdered = rbind(lgOrdered,genes)
}
lgOrdered

chrpos=0
for (i in unique(lgOrdered[,2])){
    chrpos[i+1] = max(which(lgOrdered[,2]==i))
}

dist_map = dist[is.na(dist[,7]) ==F,]
genes = as.character(dist_map[,1])
effects = as.character(dist_map[,2])

m=matrix(nrow=6300+length(unique(genes)[unique(genes) %in% rownames(lg)==F]),ncol=6300)
colnames(m) = rownames(lgOrdered)
rownames(m) = c(rownames(lgOrdered),unique(genes)[unique(genes) %in% rownames(lgOrdered)==F])
for (i in colnames(m)){
    if (i %in% rownames(loc)){
        m[i,i] = 2
    }
    if (i %in% effects){
        m[genes[which(effects %in% i)],i] =1
    }
}
m[is.na(m)] = 0
ncol(m)
nrow(m)


plotD=matrix(ncol=3,nrow=nrow(m))
for (i in 1:6300){
    loci = which(m[i,] == 1)
    lociName = names(which(m[i,] == 1))
    associations = dist[dist[,'Gene'] %in% rownames(m)[i] & dist[,'Effect'] %in% lociName,]
    top = associations[order(associations[,'PValue']),][1,]
    plotD[i,1] = i
    if (all(is.na(top))){
            print(i)
    }
    else {
        plotD[i,2] = which(colnames(m) %in% top[,'Effect'])
        plotD[i,3] = top[,'PValue']
    }
}
colnames(plotD) = c('gene','effect','P')

plotD = as.data.frame(plotD)

quartz(width=9,height=7)
ggplot(plotD,aes(plotD[,'effect'],plotD[,'gene'])) + geom_point(aes(size=as.numeric(-1*plotD[,'P']))) + scale_size(range = c(1, 3)) + scale_y_continuous(breaks=c(chrpos[chrpos<=6200])) + scale_x_continuous(breaks=c(chrpos[chrpos<=6200]))


plotL=matrix(ncol=3,nrow=nrow(m))
for (i in 1:6300){
    loci = which(m[i,] == 2)
    lociName = names(which(m[i,] == 2))
    associations = loc[rownames(loc) %in% lociName,]
    top = associations[order(associations[,'PValue']),][1,]
    plotL[i,1] = i
    if (all(is.na(top))){
        print(i)
    }
    else {
        plotL[i,2] = i
        plotL[i,3] = as.numeric(top[,'PValue'])
    }
}
colnames(plotL) = c('gene','effect','P')

plotL = as.data.frame(plotL)

quartz(width=9,height=7)
ggplot(plotL,aes(plotL[,'effect'],plotL[,'gene'])) + geom_point(aes(size=as.numeric(-1*plotL[,'P'])),alpha=0.05) + scale_size(range = c(1, 3)) + scale_y_continuous(breaks=c(chrpos[chrpos<=6200])) + scale_x_continuous(breaks=c(chrpos[chrpos<=6200]))



plotL[,'id'] = 'local'
plotD[,'id'] = 'distant'
p = rbind(plotD,plotL)
head(p)

p$bins <- cut(p$P, breaks=c(1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3), labels=c("1e-10","1e-9","1e-8","1e-7","1e-6","1e-5","1e-4"))


cols = c('#333333','#6699FF')
names(cols) = c('distant','local')

p[,'col'] = cols[p[,'id']]

quartz(width=9,height=7)
#pdf('Figure2C_revised.pdf',width=8.7,height=7)
ggplot(p[p[,'effect']<6193 & p[,'gene']<6193,],aes(effect,gene,colour=id,size=bins)) +  geom_point(aes(size=-1*as.numeric(bins)),alpha=0.5) + scale_size(range = c(1, 4))+ scale_colour_manual(breaks = p$id,values = unique(as.character(p$col))) + scale_y_continuous(breaks=c(chrpos[chrpos<=6200])) + scale_x_continuous(breaks=c(chrpos[chrpos<=6200]))
#dev.off()

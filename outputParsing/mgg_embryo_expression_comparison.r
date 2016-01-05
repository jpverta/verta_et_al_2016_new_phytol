'''
Graphs for megagametophyte vs. embryo expression comparison and expression count reproducibility.
'''

megaPath = # Path to summarized megagametophyte counts
embryoPath = # Path to summarized embryo counts
controlPath = # Path to summarized counts of control samples (run multiple times on same and different lanes 

mgg = read.table(megaPath,row.names=1)
embryo = read.table(embryoPath,row.names=1)

mgg_samples={}
for (i in colnames(mgg)){
	mgg_samples = append(mgg_samples,unlist(strsplit(i,'\\.'))[1])
}

embryo_samples = {}
for (i in colnames(embryo)){
	name = unlist(strsplit(i,'\\.'))[1]#[grepl("B",unlist(strsplit(i,'\\.')))]
	if (length(unlist(strsplit(name,''))) < 6){
		embryo_samples = append(embryo_samples,name)
	}
	else {
		embryo_samples = append(embryo_samples,unlist(strsplit(name,'_'))[2])
	}
}

mgg_samples = sub('[A-Z]','',mgg_samples)
embryo_samples = sub('[A-Z]','',embryo_samples)

# modify column names

colnames(mgg) = mgg_samples
colnames(embryo) = embryo_samples

# order samples

order={}
for (i in mgg_samples){
	order = append(order,paste(unlist(strsplit(i,'A'))[1],'B',sep=''))
}

embryo = embryo[,order]

# check

colnames(mgg)
colnames(embryo)

# rpm calculation

mgg_rpm=data.frame()
for (i in colnames(mgg)){
	mgg_rpm[1:nrow(mgg),i] = mgg[1:nrow(mgg),i]/(sum(mgg[,i])*10^(-6))
}

embryo_rpm=data.frame()
for (i in colnames(embryo)){
	embryo_rpm[1:nrow(embryo),i] = embryo[1:nrow(embryo),i]/(sum(embryo[,i])*10^(-6))
}

# graphs

plot(mgg_rpm[,1],embryo_rpm[,1],pch=20,xlim=c(0,3000),ylim=c(0,3000),main='megagametophyte vs embryo')
points(mgg_rpm[,1],mgg_rpm[,2],pch=20,col='green')
points(embryo_rpm[,1],embryo_rpm[,2],pch=20,col='blue')

pdf('counts_mgg_vs_emb.pdf')
plot(mgg_rpm[,1],embryo_rpm[,1],pch=20,xlim=c(0,3000),ylim=c(0,3000),main='megagametophyte vs embryo',cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()
pdf('counts_mgg_vs_mgg.pdf')
plot(mgg_rpm[,1],mgg_rpm[,2],pch=20,col='green',xlim=c(0,3000),ylim=c(0,3000),main='megagametophyte vs megagametophyte',cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()
pdf('counts_emb_vs_emb.pdf')
plot(embryo_rpm[,1],embryo_rpm[,2],pch=20,col='blue',xlim=c(0,3000),ylim=c(0,3000),main='embryo vs embryo',cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()

# control runs

counts = read.table(controlPath)

rpm=data.frame()
for (i in colnames(counts)){
	rpm[1:nrow(counts),i] = counts[1:nrow(counts),i]/(sum(counts[,i])*10^(-6))
}


# graphs of PRM counts

pdf('rpkm_plot.pdf')
plot(mgg_rpm[,1],embryo_rpm[,1],pch=20,xlim=c(0,3000),ylim=c(0,3000),cex.lab=1.5,cex.axis=1.5)
points(mgg_rpm[,1],mgg_rpm[,2],pch=20,col='green')
points(embryo_rpm[,1],embryo_rpm[,2],pch=20,col='blue')
points(rpm[,2],rpm[,3],pch=20,col='gray44')
points(rpm[,1],rpm[,2],pch=20,col='red')
dev.off()

pdf('counts_between_lanes.pdf')
plot(rpm[,2],rpm[,3],pch=20,col='grey44',xlim=c(0,3000),ylim=c(0,3000),main='lane A vs lane B',cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()
pdf('counts_within_lanes.pdf')
plot(rpm[,1],rpm[,2],pch=20,col='red',xlim=c(0,3000),ylim=c(0,3000),main='lane A vs lane A',cex=1.5,cex.lab=1.5,cex.axis=1.5)
dev.off()

# mean difference betwen 1N and 2N

mmean = apply(mgg_rpm,1,mean)
emean = apply(embryo_rpm,1,mean)

plot(mmean,emean,xlim=c(0,1000),ylim=c(0,1000))

library(hexbin)
hbin = hexbin(mmean[mmean<5000],emean[emean<5000],xbins=50)
plot(hbin)
hbin = hexbin(mmean[emean>0 & emean<10 & mmean>0 & mmean<10],emean[emean>0 & mmean<10 & mmean>0 & emean<10],xbins=50)
plot(hbin)

# log counts

mgg_log = log10(mgg)
embryo_log = log10(embryo)

cont_log=data.frame()
for (i in colnames(counts)){
	cont_log[1:nrow(counts),i] = log10(counts[1:nrow(counts),i])
}

mmean = apply(mgg_log,1,mean)
emean = apply(embryo_log,1,mean)
mmean[mmean == '-Inf'] = 0
emean[emean == '-Inf'] = 0
plot(mmean,emean,xlim=c(0,10),ylim=c(0,10))

fit = lm(emean~mmean)
summary(fit)

library(hexbin)
pdf('log_mean_mgg_vs_emb.pdf')
hbin = hexbin(mmean[mmean > 0 & emean> 0],emean[mmean > 0 & emean> 0],xbins=50)
plot(hbin,main='mean megagametophyte vs mean embryo')
dev.off()

pdf('log_counts_mgg_vs_emb.pdf')
hbin = hexbin(mgg_log[mgg_log[,1]>0 & embryo_log[,1]>0,1],embryo_log[mgg_log[,1]>0 & embryo_log[,1]>0,1])
plot(hbin,main='megagametophyte vs embryo')
dev.off()
pdf('log_counts_mgg_vs_mgg.pdf')
hbin = hexbin(mgg_log[mgg_log[,1]>0 & mgg_log[,2]>0,1],mgg_log[mgg_log[,1]>0 & mgg_log[,2]>0,2])
plot(hbin,main='megagametophyte vs megagametophyte')
dev.off()
pdf('log_counts_emb_vs_emb.pdf')
hbin = hexbin(embryo_log[embryo_log[,1]>0 & embryo_log[,2]>0,1],embryo_log[embryo_log[,1]>0 & embryo_log[,2]>0,2])
plot(hbin,main='emrbyo vs embryo')
dev.off()
pdf('log_counts_between_lanes.pdf')
hbin = hexbin(cont_log[cont_log[,2]>0 & cont_log[,3]>0,2],cont_log[cont_log[,2]>0 & cont_log[,3]>0,3])
plot(hbin,main='lane A vs lane B')
dev.off()
pdf('log_counts_within_lanes.pdf')
hbin = hexbin(cont_log[cont_log[,1]>0 & cont_log[,2]>0,1],cont_log[cont_log[,1]>0 & cont_log[,2]>0,2])
plot(hbin,main='lane A vs lane A')
dev.off()


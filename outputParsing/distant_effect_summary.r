'''
Non-redudant distant associations.
'''

distPath = # output from distantLRTsummary.r

library(qvalue)

dist = read.table(distPath,header=T,stringsAsFactors=F)
dist=unique(dist)
rownames(dist) = paste(dist[,1],dist[,2],sep=':')
dist = dist[is.na(dist[,6])==F,]

effects = dist[,2]
genes=dist[,1]
R=data.frame()
for (i in unique(genes)){
	dr = dist[which(genes %in% i),]
	if (nrow(dr) > length(unique(dr[,8]))){
		t=matrix(ncol=9,nrow=1)
		colnames(t) = colnames(dr)
		t[1,1:9] = NA
		dr = dr[order(as.numeric(dr[,6])),] # start with the lowest P-value associations
		for (y in 1:nrow(dr)){
			if (dr[y,8] %in% t[,8] == F){
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

nrow(R[is.na(R[,7]),]) # non-redudant, focal gene mapped
nrow(R[is.na(R[,7])==F,]) # non-redundant, focal gene non-mapped
'''
Test for distant effect hotspots
'''

distPath = # Output from distantLRTsummary.r 
lgPath = # Output from mapFiles.r (linkage blocks)

sig = read.table(distPath,stringsAsFactors=F)
lg = read.table(mapPath)

effects = sig[,2]
genes=sig[,1]
R=data.frame()
for (i in unique(genes)){
	dr = sig[which(genes %in% i),]
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

sigPos = table(R[,8])

scounts=c()
scounts[names(sigPos)] = c(sigPos)
scounts[names(table(lg[,2]))[names(table(lg[,2])) %in% names(sigPos) == F]] = 0

# choose twelve largest blocks
scounts = scounts[1:12]

# expected number of effects per linkage group:
seffects = R[R[,8] %in% c(1:12),2]
expected = length(seffects)*c(table(lg[lg[,2] %in% c(1:12),2])/nrow(lg[lg[,2] %in% c(1:12),]))

# difference between observed and expected numbers of effects
difference = scounts[names(table(lg[,2]))]-length(seffects)*c(table(lg[,2])/nrow(lg))

#test
observed = scounts
chisq.test(observed,expected)


# results:
#	Pearson's Chi-squared test
#
#data:  observed and expected
#X-squared = 120, df = 110, p-value = 0.2421

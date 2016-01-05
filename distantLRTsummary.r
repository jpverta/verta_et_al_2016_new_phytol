
'''
Summarize distant effect tests to a single table. Loops through all test outcomes generated with distantLRT.r.
'''


testDir = # PATH TO DIRECTORY WITH RESULTS FROM DISTANT EFFECT TESTS (output from distantLRT.r)

lgPath = # PATH TO FILE THAT DEFINES LINKAGE GROUPS (output from mapFiles.r)

rfPath = # PATH TO FILE THAT DESCRIBES PAIR-WISE RECOMBINAITON FRACTIONS (output from mapFiles.r)

outPath = # PATH TO OUTFILE

lg = read.table(lgPath,header=T)
rf = read.table(rfPath)

f = list.files(testDir)

results = matrix(ncol=9,nrow=length(f)*6258)
colnames(results) = c('target','effect','logFC','logCPM','LR','Pvalue','LG_target','LG_effect','RF')
results=as.data.frame(results)
for (i in f){
	trans = read.table(paste(testDir,i,sep=''),header=T)
	trans[,5:9] = NA
	trans[,3:6] = trans[,1:4]
	trans[,1] = i
	trans[,2] = rownames(trans)
	colnames(trans) = c('target','effect','logFC','logCPM','LR','Pvalue','LG_target','LG_effect','RF')
	trans[,7] = lg[i,2]
	trans[,8] = lg[trans[,2],2]
	trans[trans[,2] %in% colnames(rf),9] = as.numeric(rf[i,trans[trans[,2] %in% colnames(rf),2]])
	rownames(trans) = paste(i,rownames(trans),sep=':')
	lastrow = unlist(which(is.na(results[,1])))[1]
	results[c(lastrow):c(lastrow+nrow(trans)-1),1:9] = trans[1:nrow(trans),1:9]
}

write.table(results,outPath)

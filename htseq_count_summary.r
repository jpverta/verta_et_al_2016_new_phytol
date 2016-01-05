counts=data.frame()
for (i in list.files('.')){
	if (grepl('.htseq_counts', i)==T){
		ht = read.table(i,row.names=1)
		counts[1:nrow(ht),i] = ht[,1]
		rownames(counts) = rownames(ht)
	}
}
head(counts)

write.table(counts,'htseq_counts.txt')

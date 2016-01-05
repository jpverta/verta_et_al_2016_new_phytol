
'''
Generation of files required for haploid genetic map. Estimation of genetic map with R/qtl.
'''

snpPath = 'snp.filtered_all_combinedBlocks' # PATH TO haplotype_parser_v3.py OUTPUT
haploPath = 'snp.filtered_haplotypes_combinedBlocks' # PATH TO snp_parser_v2.py OUTPUT

'''
Lines (SNPs) from snp file and haplotype file will be parsed to select a non-redundant set of SNP genotypes. Preference will be give for genotypes as defined in haplotype file.
'''

genPath = 'genfile_FINAL'# PATH WHERE "genfile" (input to R/qtl) WILL BE SAVED

chridPath = 'chridfile_FINAL' # PATH WHERE "chridfile" (input to R/qtl) WILL BE SAVED

mnamesPath = 'mnamesfile_FINAL'# PATH WHERE "mnamesfile" (input to R/qtl) WILL BE SAVED

mapPath = # PATH WHERE FILE DESCRIBING LINKAGE GROUPS WILL BE SAVED

rfPath = # PATH WHERE FILE DESCRIBING PAIRWISE RECOMBINATION FRACTIONS WILL BE SAVED

phePath ='phefile_rpkm_all'

snp = read.table(snpPath,sep='\t',header=FALSE)
haplo = read.table(haploPath,sep='\t')
rownames(haplo) = haplo[,1]

get_samples = function(x){
         return(unlist(strsplit(unlist(strsplit(as.character(x),'],'))[2],'\\['))[2])
        }

genes = unique(c(as.character(snp[,1]),as.character(haplo[,1])))
geno = matrix(nrow=length(genes),ncol=1)
rownames(geno) = genes
for (i in rownames(geno)){
        if (i %in% rownames(haplo)){
                geno[i,1] = as.character(get_samples(haplo[i,7]))
        }
        else {
                if (i %in% snp[,1]){
                        geno[i,1] = as.character(snp[,6][snp[,1] %in% i][1])

                }
        }
}

samples = unlist(strsplit(as.character(snp[,6]),','))
samples = unique(gsub("\\]",'',gsub("\\s",'',gsub("\\[",'',samples))))
samples = unlist(strsplit(samples,'A.vcf'))

genfile = matrix(nrow=66,ncol=nrow(geno))
rownames(genfile) = samples
for (i in 1:nrow(geno)){
	for (y in samples){
		if (grepl(y,as.character(geno[i,1]))){
			genfile[y,i] = 1 }
		else {
			genfile[y,i] = 0 }
	}		
}
rownames(genfile) = gsub("[a-z_.,]",'',rownames(genfile))
colnames(genfile) = rownames(geno)
write.table(genfile,genPath,row.names=FALSE,col.names=FALSE)

chridfile = rep(1, times=ncol(genfile))
mnamesfile = paste(rownames(geno),sep=';')
write.table(chridfile,chridPath,row.names=FALSE,col.names=FALSE)
write.table(mnamesfile,mnamesPath,row.names=FALSE,col.names=FALSE)

# GENETIC MAP

library(qtl)

cross = read.cross(format=c("gary"),estimate.map=FALSE,genotypes=c("AA","BB"),alleles=c("R","A"),genfile=genPath,mapfile=NULL,phefile=phePath,chridfile=chridPath,mnamesfile=mnamesPath)
class(cross)[1] <- "dh"

cross = formLinkageGroups(cross, max.rf=0.35, min.lod=6,reorgMarkers=TRUE)
plot(cross)
toswitch = markernames(cross, chr=c(4,8,11,13,14,16,17,18,21,22,23,24,25))
cross = switchAlleles(cross, toswitch)

cross2 = est.rf(cross2)
cross2 = formLinkageGroups(cross2, max.rf=0.35, min.lod=6,reorgMarkers=TRUE)
plotRF(cross2)

map = pull.map(cross2)
LG=data.frame()
for (lg in names(map)){
	for (gene in names(map[[lg]])){
		LG[gene,1] = 1
		LG[gene,2] = lg
	}
}	

rf = pull.rf(cross2,what='rf')

write.table(map,lgPath)
write.table(rf,rfPath)

# Order markers on linkage groups
cross = orderMarkers(cross)


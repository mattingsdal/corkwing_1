setwd("C:/Users/mortenma/Google Drive/post doc/essential files/paper 1/test")

library(SNPRelate)
library(pheatmap)
library(Cairo)

bed.fn <- "SNP.QC.population.bed"
fam.fn <- "SNP.QC.population.fam"
bim.fn <- "SNP.QC.population.bim"
snpgdsBED2GDS(bed.fn, fam.fn, bim.fn, "SNP.QC.population.gds")
snpgdsSummary("SNP.QC.population.gds")

genofile <- snpgdsOpen("SNP.QC.population.gds")
id <- read.gdsn(index.gdsn(genofile, "sample.id"))

## Fst, F
pop=c(rep("F",7),rep("A",8),rep("E",8),rep("H",10),rep("C",8),rep("B",8),rep("D",8),rep("G",8))
pop2=c(rep("AR",7),rep("UK",8),rep("EG",8),rep("GM",10),rep("NH",8),rep("SM",8),rep("ST",8),rep("TV",8))

res=matrix(ncol=8,nrow=8)
index=unique(pop)

for (j in 1:8){
for (i in 1:8){
	subpop=pop[pop%in%index[j]|pop%in%index[i]]
	subsample=id[pop%in%index[j]|pop%in%index[i]]
	res[i,j]=ifelse(length(unique(subpop))==2,snpgdsFst(genofile, sample.id=subsample, population=as.factor(subpop),method="W&C84",autosome.only=FALSE,verbose=FALSE)$Fst,c(0))
	}
}
colnames(res)= unique(pop2)
row.names(res)= unique(pop2)

FST=res

library(pheatmap)
library(RColorBrewer)

Cairo(file="FST_heatmap.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = rep(0, 4))
pheatmap(FST,border="white",color = brewer.pal(9,"Blues"),cluster_rows=TRUE,cluster_cols=TRUE,display_numbers=T,number_color="gray60",number_format="%.3f",fontsize_number=3,fontsize_row=6,fontsize_col=6,legend=FALSE,treeheight_row=25,treeheight_col=25)
dev.off()



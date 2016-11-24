plink --bfile SNP.QC.population --freq --keep plink/pop/south --allow-extra-chr --out south
plink --bfile SNP.QC.population --freq --keep plink/pop/west --allow-extra-chr --out west
plink --bfile SNP.QC.population --freq --keep plink/pop/ard --allow-extra-chr --out ard

R
library(limma)
library(Cairo)

setwd("C:/Users/mortenma/Google Drive/post doc/essential files/paper 1/Figures")
south=read.table("plink_south_imputed.frq",sep="",header=T)
west=read.table("plink_west_imputed.frq",sep="",header=T)
ard=read.table("plink_ard_imputed.frq",sep="",header=T)

venn=matrix(ncol=3,nrow=nrow(south))
venn[,1]=ifelse(as.numeric(west[,5])>0,1,0)
venn[,2]=ifelse(as.numeric(south[,5])>0,1,0)
venn[,3]=ifelse(as.numeric(ard[,5])>0,1,0)

colnames(venn)=c("West","South","Ard")
row.names(venn)=south[,2]

venn2=venn[!rowSums(venn)==0,]

a <- vennCounts(venn2)
Cairo(file="vennDiagram.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = c(0, 0,0,0))
vennDiagram(a,circle.col=c("#ff4500","#4169e1","#ffd700"),cex=0.5,lwd=1,mar=c(0,0,0,0))
dev.off()


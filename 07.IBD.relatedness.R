setwd("C:/Users/mortenma/Google Drive/post doc/essential files/paper 1/test")

library(SNPRelate)
library(pheatmap)
library(Cairo)

genofile <- snpgdsOpen("SNP.QC.population.gds")
id <- read.gdsn(index.gdsn(genofile, "sample.id"))
snpid <- read.gdsn(index.gdsn(genofile, "snp.id"))

snpid=
pop=c(rep("F",7),rep("A",8),rep("E",8),rep("H",10),rep("C",8),rep("B",8),rep("D",8),rep("G",8))
pop2=c(rep("AR",7),rep("UK",8),rep("EG",8),rep("GM",10),rep("NH",8),rep("SM",8),rep("ST",8),rep("TV",8))

ibd <- snpgdsIBDMLE(genofile, sample.id=id,snp.id=snpid,maf=0.05, missing.rate=0.05, num.thread=2,autosome.only=FALSE)

pop2=c(rep("AR",7),rep("UK",8),rep("EG",8),rep("GM",10),rep("NH",8),rep("SM",8),rep("ST",8),rep("TV",8))

ibd2=snpgdsIBDSelection(ibd)

Cairo(file="IBD_relatedness.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = rep(2, 4))
plot(ibd2$k0,ibd2$k1,ylim=c(0,1),xlim=c(0,1),xlab="k0",ylab="k1",cex.axis=0.5,lwd=0.5,pch=19,cex=0.1)
lines(c(0,1), c(1,0), col="black", lty=2)
ibdAreasDraw(xcol=rep("black",6),rel.lwd=0.5)
text(0.25,0.5,label=c("1st degree"),col="black",cex=0.3)
text(0.6,0.5,label=c("2nd degree"),col="black",cex=0.3)
text(0.8,0.3,label=c("3rd degree"),col="black",cex=0.3)
text(0.93,0.15,label=c("unrelated"),col="black",cex=0.3)
text(0.16,0.95,label=c("duplicate"),col="black",cex=0.3)
text(0.12,0.71,"ST46 & ST48",cex=0.2)
text(0.45,0.45,"all UK samples",cex=0.2)
dev.off()


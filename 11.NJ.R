library("adegenet")
library(ape)
library(Cairo)

setwd("C:/Users/mortenma/Google Drive/post doc/essential files/population structure")

data=read.PLINK("SNP.QC.population.raw")
pop(data)=c(rep(3,7),rep(1,8),rep(3,18),rep(2,24),rep(3,8))


## NJ tree of distance
tre <- nj(dist(as.matrix(data)))


pca1 <- glPca(data)

## nice PCA plot
#s.class(pca1$scores, pop(res),col=c("#ffd700","#ff4500","#4169e1"),clabel="NULL")

# default PCA plot
tmp=c(rep("F",7),rep("A",8),rep("E",8),rep("H",10),rep("C",8),rep("B",8),rep("D",8),rep("G",8))
tmp2=c(rep("#4169e1",7),rep("#ffd700",8),rep("#4169e1",8),rep("#4169e1",10),rep("#ff4500",8),rep("#ff4500",8),rep("#ff4500",8),rep("#4169e1",8))
pca2=cbind(pca1$scores,tmp,tmp2)

# plot
Cairo(file="NJ_UNROOTED.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = rep(0, 4))
plot(tre, typ="unrooted", cex=0.2,label.offset=5,no.margin=TRUE,lab4ut="axial",font=2)
tiplabels(pch=20, col=scales::alpha(pca2[,4],0.5), cex=1)
dev.off()

Cairo(file="NJ_FAN.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = rep(0, 4))
plot(tre, typ="fan", cex=0.2,label.offset=5,no.margin=TRUE,lab4ut="axial",font=2)
tiplabels(pch=20, col=scales::alpha(pca2[,4],0.5), cex=1)
dev.off()


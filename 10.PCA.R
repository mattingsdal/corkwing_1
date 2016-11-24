library(adegenet)
library(ape)
library(Cairo)
library(car)

setwd("C:/Users/mortenma/Google Drive/post doc/essential files/paper 1/test")

data=read.PLINK("SNP.QC.population.raw",parallel=FALSE)
pop(data)=c(rep(3,7),rep(1,8),rep(3,18),rep(2,24),rep(3,8))

pca1 <- glPca(data,parallel=FALSE,nf=2)


# all samples
tmp=c(rep("F",7),rep("A",8),rep("E",8),rep("H",10),rep("C",8),rep("B",8),rep("D",8),rep("G",8))
tmp2=c(rep("#4169e1",7),rep("#ffd700",8),rep("#4169e1",8),rep("#4169e1",10),rep("#ff4500",8),rep("#ff4500",8),rep("#ff4500",8),rep("#4169e1",8))
plotpca1=cbind(pca1$scores,tmp,tmp2)

Cairo(file="PCA.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = c(1, 1,1,1))
plot(plotpca1[,1:2],pch=19,cex=1,cex.axis=0.3,cex.lab=0.3,lwd=0,xaxt='n',yaxt='n',col=scales::alpha(plotpca1[,4],0.5))
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
text(plotpca1[,1:2],label=as.character(plotpca1[,3]),cex=0.3,font=2,col=scales::alpha("black",1))
limit1=quantile(range(as.numeric(plotpca1[,1])))
limit2=quantile(range(as.numeric(plotpca1[,2])))
pc1=pca1$eig[1]/sum(pca1$eig)*100
pc2=pca1$eig[2]/sum(pca1$eig)*100
text(limit1[2],limit2[1],label=paste("PC1=",round(pc1,2),"%",sep=""),cex=0.5)
text(limit1[4],limit2[1],label=paste("PC2=",round(pc2,2),"%",sep=""),cex=0.5)
dev.off()

## exclude scotland
data2=data[-8:-15,]
pca2 <- glPca(data2,parallel=FALSE,nf=2)
tmp=c(rep("F",7),rep("E",8),rep("H",10),rep("C",8),rep("B",8),rep("D",8),rep("G",8))
tmp2=c(rep("#4169e1",7),rep("#4169e1",8),rep("#4169e1",10),rep("#ff4500",8),rep("#ff4500",8),rep("#ff4500",8),rep("#4169e1",8))
plotpca2=cbind(pca2$scores,tmp,tmp2)
### only norway
Cairo(file="PCA_norway.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = c(1, 1,1,1))
plot(plotpca2[,1:2],pch=19,cex=1,cex.axis=0.3,cex.lab=0.3,lwd=0,xaxt='n',yaxt='n',col=scales::alpha(plotpca2[,4],0.5))
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
pc1=pca2$eig[1]/sum(pca2$eig)*100
pc2=pca2$eig[2]/sum(pca2$eig)*100
text(plotpca2[,1:2],label=as.character(plotpca2[,3]),cex=0.3,font=2,col=scales::alpha("black",1))
limit1=quantile(range(as.numeric(plotpca2[,1])))
limit2=quantile(range(as.numeric(plotpca2[,2])))
text(limit1[2],limit2[1],label=paste("PC1=",round(pc1,2),"%",sep=""),cex=0.5)
text(limit1[3],limit2[1],label=paste("PC2=",round(pc2,2),"%",sep=""),cex=0.5)
#dataEllipse(pca3$scores,groups=as.factor(plotpca3[,3]),robust=T,plot.points=FALSE,add=TRUE,col=c("#4169e1","#6941e1","#41b9e1","#839deb"),lwd=0.5,center.pch=F,ellipse.label=F,levels=c(0.8))
dev.off()

## only south
south=data2[-26:-49,]
south$pop=as.factor(c(rep(1,7),rep(2,8),rep(3,10),rep(4,8)))
pca3 <- glPca(south,parallel=FALSE,nf=2)
tmp=c(rep("F",7),rep("E",8),rep("H",10),rep("G",8))
tmp2=c(rep("#4169e1",7),rep("#6941e1",8),rep("#41b9e1",8),rep("#839deb",10))
plotpca3=cbind(pca3$scores,tmp,tmp2)
Cairo(file="PCA_south.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = c(1, 1,1,1))
plot(plotpca3[,1:2],pch=19,cex=1,cex.axis=0.3,cex.lab=0.3,lwd=0,xaxt='n',yaxt='n',col=scales::alpha(plotpca3[,4],0.5))
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
pc1=pca3$eig[1]/sum(pca3$eig)*100
pc2=pca3$eig[2]/sum(pca3$eig)*100
text(plotpca3[,1:2],label=as.character(plotpca3[,3]),cex=0.3,font=2,col=scales::alpha("black",1))
limit1=quantile(range(as.numeric(plotpca3[,1])))
limit2=quantile(range(as.numeric(plotpca3[,2])))
text(limit1[2],limit2[1],label=paste("PC1=",round(pc1,2),"%",sep=""),cex=0.5)
text(limit1[3],limit2[1],label=paste("PC2=",round(pc2,2),"%",sep=""),cex=0.5)
#dataEllipse(pca3$scores,groups=as.factor(plotpca3[,3]),robust=T,plot.points=FALSE,add=TRUE,col=c("#4169e1","#6941e1","#41b9e1","#839deb"),lwd=0.5,center.pch=F,ellipse.label=F,levels=c(0.8))
dev.off()

## only west
west=data2[26:49,]
west$pop=as.factor(c(rep(1,8),rep(2,8),rep(3,8)))
pca4 <- glPca(west,parallel=FALSE,nf=2)
tmp=c(rep("B",8),rep("C",8),rep("D",8))
tmp2=c(rep("#b33000",8),rep("#ff4500",8),rep("#ff7d4d",8))
plotpca4=cbind(pca4$scores,tmp,tmp2)
Cairo(file="PCA_west.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = c(1, 1,1,1))
plot(plotpca4[,1:2],pch=19,cex=1,cex.axis=0.3,cex.lab=0.3,lwd=0,xaxt='n',yaxt='n',col=scales::alpha(plotpca4[,4],0.5))
abline(h=0,col="lightgray")
abline(v=0,col="lightgray")
pc1=pca4$eig[1]/sum(pca4$eig)*100
pc2=pca4$eig[2]/sum(pca4$eig)*100
text(plotpca4[,1:2],label=as.character(plotpca4[,3]),cex=0.3,font=2,col=scales::alpha("black",1))
limit1=quantile(range(as.numeric(plotpca4[,1])))
limit2=quantile(range(as.numeric(plotpca4[,2])))
text(limit1[2],limit2[1],label=paste("PC1=",round(pc1,2),"%",sep=""),cex=0.5)
text(limit1[3],limit2[1],label=paste("PC2=",round(pc2,2),"%",sep=""),cex=0.5)
#dataEllipse(pca4$scores,levels=.8,groups=as.factor(plotpca4[,3]),robust=T,plot.points=FALSE,add=TRUE,col=c("#b33000","#ff4500","#ff7d4d"),lwd=0.5,center.pch=F,ellipse.label=F)
dev.off()





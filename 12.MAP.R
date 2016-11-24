setwd("E:/data/bam_fixrg_dedup/freebayes/popgen")

#########################################################################
#########################################################################
#### construct map of interest, Northsea
library(scales)
library(rworldmap)
library(igraph)

newmap <- getMap(resolution = "high")
plot(newmap, xlim = c(-6, 15), ylim = c(55, 63), asp = 2)
loc=matrix(ncol=2,nrow=7)

Ardtoe=c(-5.88232,56.76745)
Smola=c(8.00687,63.38657)
Norheimsund=c(6.14564,60.37089)
Stavanger=c(5.73311,58.96998)
Egersund=c(5.99980,58.45142)
Arendal=c(8.77245,58.46176)
Tvedestrand=c(8.93140,58.62228)
Gullmarn=c(11.43333,58.25)


index=c("A","B","C","D","E","F","G","H")
col=c("#ffd700","#ff4500","#ff4500","#ff4500","#4169e1","#4169e1","#4169e1","#4169e1")


tmp=data.frame(rbind(Ardtoe,Smola,Norheimsund,Stavanger,Egersund,Arendal,Tvedestrand,Gullmarn))
loc=cbind(tmp,index,col)

Cairo(file="MAP.pdf",type="pdf",width=80,height=80,units="mm")
par(mar = rep(0, 4))
plot(newmap, xlim = c(-6, 15), ylim = c(55, 63), asp = 2,lwd=1,col="gray95")
points(loc, col = alpha(loc$col,0.5), cex = 2,lwd=0,pch=19)
points(loc, col = "black", cex = 0.5,lwd=0,pch=19)

rect(loc[1,1]+0.5,loc[1,2]-0.25,loc[1,1]+1.5,loc[1,2]+0.25,col="white",border=NA)
rect(loc[2,1]+0.5,loc[2,2]-0.25,loc[2,1]+1.5,loc[2,2]+0.25,col="white",border=NA)
rect(loc[3,1]+0.5,loc[3,2]-0.25,loc[3,1]+1.5,loc[3,2]+0.25,col="white",border=NA)
rect(loc[4,1]+0.5,loc[4,2]-0.25,loc[4,1]+1.5,loc[4,2]+0.25,col="white",border=NA)
rect(loc[5,1]+0.5,loc[5,2]-0.25,loc[5,1]+1.5,loc[5,2]+0.25,col="white",border=NA)
rect(loc[6,1]+0.5,loc[6,2]-0.25,loc[6,1]+1.5,loc[6,2]+0.25,col="white",border=NA)
rect(loc[7,1]+0.5,loc[7,2]-0.25,loc[7,1]+1.5,loc[7,2]+0.25,col="white",border=NA)
rect(loc[8,1]+0.5,loc[8,2]-0.25,loc[8,1]+1.5,loc[8,2]+0.25,col="white",border=NA)

#text(x=loc[,1]+1,y=loc[,2],label=index,cex=1,font=2,col=alpha("black",1))
loc[6,2]=58.17000

text(x=loc[,1]+1,y=loc[,2],label=index,cex=1,font=2,col=alpha("black",1))

dev.off()

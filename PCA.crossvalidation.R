setwd("C:/Users/mortenma/Google Drive/post doc/essential files/paper 1/test")

library(adegenet)
library(caret)
library(doParallel)

data=read.PLINK("SNP.QC.population.raw",parallel=FALSE)

# norway only
norway=data[-8:-15,]
# only south
south=norway[-26:-49,]
# only west
west=norway[26:49,]

pop(norway)=c(rep("south",7),rep("south",8),rep("south",10),rep("west",24),rep("south",8))
pop(south)=c(rep("AR",7),rep("EG",8),rep("GM",10),rep("TV",8))
pop(west)=c(rep("NH",8),rep("SM",8),rep("ST",8))

norway.pca <- glPca(norway,nf=2,parallel=FALSE)
south.pca <- glPca(south,nf=5,parallel=FALSE)
west.pca <- glPca(west,nf=3,parallel=FALSE)

norway.load1=abs((norway.pca$loadings[,1])*(norway.pca$eig/(sum(norway.pca$eig)))[1])
norway.load2=abs((norway.pca$loadings[,2])*(norway.pca$eig/(sum(norway.pca$eig)))[2])
norway.load=norway.load1+norway.load2
names(norway.load)=norway$loc.names
norway.load=norway.load[order(norway.load,decreasing=T)]

west.load1=abs((west.pca$loadings[,1])*(west.pca$eig/(sum(west.pca$eig)))[1])
west.load2=abs((west.pca$loadings[,2])*(west.pca$eig/(sum(west.pca$eig)))[2])
west.load3=abs((west.pca$loadings[,3])*(west.pca$eig/(sum(west.pca$eig)))[3])
west.load=west.load1+west.load2+west.load3
names(west.load)=west$loc.names
west.load=west.load[order(west.load,decreasing=T)]

south.load1=abs((south.pca$loadings[,1])*(south.pca$eig/(sum(south.pca$eig)))[1])
south.load2=abs((south.pca$loadings[,2])*(south.pca$eig/(sum(south.pca$eig)))[2])
south.load3=abs((south.pca$loadings[,3])*(south.pca$eig/(sum(south.pca$eig)))[3])
south.load4=abs((south.pca$loadings[,4])*(south.pca$eig/(sum(south.pca$eig)))[4])
south.load5=abs((south.pca$loadings[,5])*(south.pca$eig/(sum(south.pca$eig)))[5])
south.load=south.load1+south.load2+south.load3+south.load4+south.load5
names(south.load)=south$loc.names
south.load=south.load[order(south.load,decreasing=T)]

### now we have summed up loadings for the most important compnenets from each of the PCA's
### now run through the top SNPs and test if classification methods can pretict status

# prepare output
result.norway.knn=matrix(ncol=7,nrow=100)
result.south.knn=matrix(ncol=7,nrow=100)
result.west.knn=matrix(ncol=7,nrow=100)
result.norway.2.knn=matrix(ncol=7,nrow=100)

en=c(1:9)
to=seq(from = 10, to = 100, by = 10)
tre=seq(from = 200, to = 1000, by = 100)

myindex=c(en,to,tre)


# define cross validation scheme
train_control <- trainControl(method="repeatedcv", number=5, repeats=10,classProbs=TRUE)

cl <- makeCluster(4)
registerDoParallel(cl)

	for (j in 1:100){
	
	## south vs west
	test.norway=cbind(norway$pop,as.data.frame(norway[,norway$loc.names%in%names(norway.load[1:j])]))
	colnames(test.norway)[1]="pop"
	
	test.south=cbind(south$pop,as.data.frame(south[,south$loc.names%in%names(south.load[1:j])]))
	colnames(test.south)[1]="pop"

	test.west=cbind(west$pop,as.data.frame(west[,west$loc.names%in%names(west.load[1:j])]))
	colnames(test.west)[1]="pop"
	
	# all sub-regions
	test.norway.2=test.norway
	test.norway.2$pop=c(rep("AR",7),rep("EG",8),rep("GM",10),rep("NH",8),rep("SM",8),rep("ST",8),rep("TV",8))
	

	# KNN
	model.norway.knn <- train(as.factor(pop) ~ ., data=test.norway, trControl=train_control, method="kknn")
	model.norway.2.knn <- train(as.factor(pop) ~ ., data=test.norway.2, trControl=train_control, method="kknn")
	model.south.knn <- train(as.factor(pop) ~ ., data=test.south, trControl=train_control, method="kknn")
	model.west.knn <- train(as.factor(pop) ~ ., data=test.west, trControl=train_control, method="kknn")
	
	knn3Train 
	result.norway.knn[j,]=rbind(as.matrix(model.norway.knn$result[1,]))
	result.norway.2.knn[j,]=rbind(as.matrix(model.norway.2.knn$result[1,]))
	result.west.knn[j,]=rbind(as.matrix(model.west.knn$result[1,]))
	result.south.knn[j,]=rbind(as.matrix(model.south.knn$result[1,]))
	
}	
stopCluster(cl); registerDoSEQ();

colnames(result.norway.knn)=names(model.norway.knn$result)
colnames(result.norway.2.knn)=names(model.norway.2.knn$result)
colnames(result.west.knn)=names(model.west.knn$result)
colnames(result.south.knn)=names(model.south.knn$result)

library(Hmisc)

### make plots
x=1:100
## KNN
plot(result.norway.knn[,4],ylim=c(0,1),type="b",pch=15,cex=1,xlab="# SNPs",ylab="Accuracy with SD",main="KNN")
arrows(x, as.numeric(result.norway.knn[,4])-as.numeric(result.norway.knn[,6]), x, as.numeric(result.norway.knn[,4])+as.numeric(result.norway.knn[,6]), length=0.05, angle=90, code=0)
points(result.norway.2.knn[,4],type="b",col="blue",pch=15,cex=1)
arrows(x, as.numeric(result.norway.2.knn[,4])-as.numeric(result.norway.2.knn[,6]), x, as.numeric(result.norway.2.knn[,4])+as.numeric(result.norway.2.knn[,6]), length=0.05, angle=90, code=0,col="blue")
points(result.west.knn[,4],type="b",col="red",pch=15,cex=1)
arrows(x, as.numeric(result.west.knn[,4])-as.numeric(result.west.knn[,6]), x, as.numeric(result.west.knn[,4])+as.numeric(result.west.knn[,6]), length=0.05, angle=90, code=0,col="red")
points(result.south.knn[,4],type="b",col="green",pch=15,cex=1)
arrows(x, as.numeric(result.south.knn[,4])-as.numeric(result.south.knn[,6]), x, as.numeric(result.south.knn[,4])+as.numeric(result.south.knn[,6]), length=0.05, angle=90, code=0,col="green")
rect(-10,1.01,150,1.2,col="white")



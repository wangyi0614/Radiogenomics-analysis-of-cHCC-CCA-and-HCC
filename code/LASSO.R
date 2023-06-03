#LASSO
library("glmnet")
library("survival")

setwd("C:\\Users\\Lenovo\\Desktop\\chccfigure")                
rt=read.table("tpm.txt",header=T,sep="\t",row.names=1)          
rt$futime=rt$futime/12

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$futime,rt$fustat))
fit=glmnet(x, y, family = "cox", maxit = 1000)
cvfit=cv.glmnet(x, y, family="cox", maxit = 1000)

coef=coef(fit, s = cvfit$lambda.min)
index=which(coef != 0)
actCoef=coef[index]
lassoGene=row.names(coef)[index]
geneCoef=cbind(Gene=lassoGene,Coef=actCoef)
write.table(geneCoef,file="geneCoef.txt",sep="\t",quote=F,row.names=F)

trainFinalGeneExp=rt[,lassoGene]
myFun=function(x){crossprod(as.numeric(x),actCoef)}
trainScore=apply(trainFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(trainScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(trainScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="TrainRisk.txt",sep="\t",quote=F,row.names=F)

rt=read.table("chcc2.txt",header=T,sep="\t",row.names=1)
rt[,3:ncol(rt)][rt[,3:ncol(rt)]<0]=0
rt$futime=rt$futime/12
testFinalGeneExp=rt[,lassoGene]
testScore=apply(testFinalGeneExp,1,myFun)
outCol=c("futime","fustat",lassoGene)
risk=as.vector(ifelse(testScore>median(trainScore),"high","low"))
outTab=cbind(rt[,outCol],riskScore=as.vector(testScore),risk)
write.table(cbind(id=rownames(outTab),outTab),file="validationRisk.txt",sep="\t",quote=F,row.names=F)

library(ggplot2)
library(gplots)
library(gridExtra)
library(car)
library(vegan)
library(reshape)
library(reshape2)
library(broom)
library(pvclust)
library(dendextend)
library(colorspace)
library(circlize)

rm(list=ls())

#######   Read in data from master folder holding all sample folders - - - - - - - - - - - - - - - -
#################
#################

col<-c(0,1,2,3)
stage<-c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28)
metric<-c("Quant", "Sigma", "Connect", "Tau1", "Tau2", "NPQ", "qP", "ABQ","qm")
metrix<-c(1,2,3,4,5,6,7,8,9)
h<-1
metaphys<-list()
for( I in 1:4)
{
	for (TT in 1:28)
	{
		for (FF in 1:9)
		{
			uni<-paste(metric[FF], stage[TT], col[I], sep=".")
			metaphys[[h]]<-data.frame("uni"=uni, "metric"=metric[FF], "metrix"=metrix[FF], "stage"=stage[TT], "col"=(col[I]+1))
			h<-h+1
		}
	}
}
metafluor1<-data.frame(do.call(rbind, metaphys))
dim(metafluor1)

#################
#################
#################

FullPheno<-read.csv("FullPheno.csv", header=TRUE)
Aller<-read.csv("DataMoteFull.csv", header=TRUE)
rownames(Aller)<-Aller[,1]
Alls<-Aller[,-1]

rownames(FullPheno)<-FullPheno[,1]
FullPheno<-FullPheno[,-1]
dim(FullPheno)

##############
##############
#####    #####
#####    #####
#####    #####
#####    #####
#####    #####
##############
##############

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

samp<-paste("MOTE-BleachingResponse", ".pdf", sep="")
pdf(file = samp, width = 6, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(2, 2.5, 0.1, 3.1), mar = c(1.1, 0.1, 1.1, 2.1))

layoutmatrix<-matrix(c(1,2,3,3), ncol=2)
layout(layoutmatrix)

hist(Allstart[,"ChaQuant"], breaks=40, axes=F, main="", xlim=c(-100, 150), ylim=c(0,15))
axis(2, at=, col.axis="black", las=2, cex.axis=1.25)
box(col="black")
hist(Allstart[,"Cha675"], breaks=30, axes=F, main="", xlim=c(-100, 150), ylim=c(0,15))
axis(1, at=, col.axis="black", las=1, cex.axis=1.25)
axis(2, at=, col.axis="black", las=2, cex.axis=1.25)
box(col="black")

plot(Allstart[,"ChaQuant"], Allstart[,"Cha675"], pch=19, xlim=c(-150,150),ylim=c(-150,150), axes=F)
axis(1, at=, col.axis="black", las=1, cex.axis=1.25)
axis(4, at=, col.axis="black", las=2, cex.axis=1.25)
box(col="black")
abline(lm(Allstart[,"Cha675"] ~ Allstart[,"ChaQuant"]), col="red", lwd=1.5)
summary(lm(Allstart[,"Cha675"] ~ Allstart[,"ChaQuant"]))
cor.test(Allstart[,"Cha675"], Allstart[,"ChaQuant"], method=c("spearman"))

dev.off()


write.csv(rownames(Allstart), "MOTEnames.csv")

shapiro.test(Allstart[,"ChaQuant"])
mean(Allstart[,"ChaQuant"])
sd(Allstart[,"ChaQuant"])

shapiro.test(Allstart[,"Cha675"])
mean(Allstart[,"Cha675"])
sd(Allstart[,"Cha675"])

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

library(caret)

Alls.cntn<-Allstart[,1:925]
Alls.trans<-preProcess(x=Alls.cntn,method=c("center","scale"))
Alls.trans

Alls.preproc<-predict(Alls.trans,newdata=Alls.cntn)
Alls.preproc<-subset(Alls.preproc, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
#Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.2, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.2, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]

dim(Alls.preproc)

descrCor<-cor(Alls.preproc)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.85)
filteredDescr<-(Alls.preproc[,-highlyCorDescri])
colnames(filteredDescr)
dim(filteredDescr)

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

Asllls<-data.frame(filteredDescr, "ChaQuant"=Allstart[,"ChaQuant"], "Cha675"=Allstart[,"Cha675"])
dim(Asllls)
lmcor675<-list()
lmcorQuant<-list()
rhocor675<-list()
rhocorQuant<-list()
for (R in 1:ncol(filteredDescr))
{
	lmcor675[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=summary(lm(filteredDescr[,R] ~ Allstart[,"Cha675"]))[[9]], "pval"=lmp(lm(filteredDescr[,R] ~ Allstart[,"Cha675"])))
	lmcorQuant[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=summary(lm(filteredDescr[,R] ~ Allstart[,"ChaQuant"]))[[9]], "pval"=lmp(lm(filteredDescr[,R] ~ Allstart[,"ChaQuant"])))
	rhocor675[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=cor.test(filteredDescr[,R], Allstart[,"Cha675"], method=c("pearson"))[4][[1]], "pval"=cor.test(filteredDescr[,R], Allstart[,"Cha675"], method=c("pearson"))[3][[1]])
	rhocorQuant[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=cor.test(filteredDescr[,R], Allstart[,"ChaQuant"], method=c("pearson"))[4][[1]], "pval"=cor.test(filteredDescr[,R], Allstart[,"ChaQuant"], method=c("pearson"))[3][[1]])
}
#lm675<-data.frame(do.call(rbind, lmcor675))
lm675<-data.frame(do.call(rbind, rhocor675))
#lmQuant<-data.frame(do.call(rbind, lmcorQuant))
lmQuant<-data.frame(do.call(rbind, rhocorQuant))

result675<-data.frame(metafluor1[match(lm675$metrics, metafluor1$uni),], lm675)
result675<-subset(result675, metric != "Connect")
result675<-result675[order(abs(result675$r2), decreasing=TRUE),]
result675<-subset(result675, pval<0.05)
dim(result675)

resultQuant<-data.frame(metafluor1[match(lmQuant$metrics, metafluor1$uni),], lmQuant)
resultQuant<-subset(resultQuant, metric != "Connect")
resultQuant<-resultQuant[order(abs(resultQuant$r2), decreasing=TRUE),]
resultQuant<-subset(resultQuant, pval<0.05)
dim(resultQuant)

write.csv(result675, "abs675.csv")
write.csv(resultQuant, "chaQuant.csv")
result675<-read.csv("abs675.csv", header=TRUE)
resultQuant<-read.csv("chaQuant.csv", header=TRUE)
dim(result675)
dim(resultQuant)

light<-c(0,0,200,200,200,200,200,200,200,200,94,94,94,94,94,94,400,400,400,400,400,400,0,0,0,0,0,0)

samp<-paste("MOTE-correlationData", ".pdf", sep="")
pdf(file = samp, width = 5, height = 3, bg="transparent")

par(mfrow=c(2,2), oma = c(2, 2.5, 0.1, 0.1), mar = c(0.1, 0.1, 0.1, 1.5))

caa<-resultQuant
ca1Rich<-data.frame(metafluor1[match(caa$metrics, metafluor1$uni),], caa)
ca1Rich1<-data.frame(aggregate(ca1Rich[,2], by=list(ca1Rich$metric), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$metric), FUN=mean))
ca1m<-data.frame(ca1Rich1, "percent"=ca1Rich1$x/dim(ca1Rich)[1])
ca1Rich2<-data.frame(aggregate(ca1Rich[,5], by=list(ca1Rich$col), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$col), FUN=mean))
ca1c<-data.frame(ca1Rich2, "percent"=ca1Rich2$x/dim(ca1Rich)[1])
ca1Rich3<-data.frame(aggregate(ca1Rich[,4], by=list(ca1Rich$stage), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$stage), FUN=mean))
ca1t<-data.frame(ca1Rich3, "percent"=ca1Rich3$x/dim(ca1Rich)[1])

tttQ<-ca1t[,1:2]

plot(light, type="l", ylim=c(0,3000), axes=F, col="gray", col.axis="gray")
axis(4, at=, col.axis="black", las=2)
box(col="black")
par(new=TRUE)
plot(tttQ$Group.1, tttQ$x, type="l", axes=F, ylim=c(-10,30))
axis(2, at=, col.axis="black", las=2)
box(col="black")

slices<-ca1m$x
lbls<-ca1m$Group.1
pie(slices, labels=lbls, col=c("#01665e", "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a")[as.factor(lbls)], cex=0.00001)

caa<-result675
ca1Rich<-data.frame(metafluor1[match(caa$metrics, metafluor1$uni),], caa)
ca1Rich1<-data.frame(aggregate(ca1Rich[,2], by=list(ca1Rich$metric), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$metric), FUN=mean))
ca1m<-data.frame(ca1Rich1, "percent"=ca1Rich1$x/dim(ca1Rich)[1])
ca1Rich2<-data.frame(aggregate(ca1Rich[,5], by=list(ca1Rich$col), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$col), FUN=mean))
ca1c<-data.frame(ca1Rich2, "percent"=ca1Rich2$x/dim(ca1Rich)[1])
ca1Rich3<-data.frame(aggregate(ca1Rich[,4], by=list(ca1Rich$stage), FUN=NROW), aggregate(abs(ca1Rich[,13]), by=list(ca1Rich$stage), FUN=mean))
ca1t<-data.frame(ca1Rich3, "percent"=ca1Rich3$x/dim(ca1Rich)[1])

ttt675<-ca1t[,1:2]
samp675<-c(7,11,12,27)
t675<-c(0,0,0,0)
QTT675<-rbind(ttt675, cbind("Group.1"=samp675, "x"=t675))
QTT675<-QTT675[order(QTT675$Group.1),]

plot(light, type="l", ylim=c(0,3000), axes=F, col="gray", col.axis="gray")
axis(4, at=, col.axis="black", las=2)
box(col="black")
par(new=TRUE)
plot(QTT675$Group.1, QTT675$x, type="l", axes=F, ylim=c(-10,30))
axis(2, at=, col.axis="black", las=2)
axis(1, at=, col.axis="black", las=1)
box(col="black")

slices<-ca1m$x
lbls<-ca1m$Group.1
pie(slices, labels=lbls, col=c("#01665e", "#35978f", "#80cdc1", "#c7eae5", "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a")[as.factor(lbls)], cex=0.00005)

dev.off()

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
library(caret)
library(randomForest)
library(Boruta)

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

bleachm<-"Cha675"
Alls.cntn<-Allstart[,1:916]
Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]
dim(Alls.preproc)
Alls.trans<-preProcess(x=Alls.preproc,method=c("center","scale"))
Alls.preproc2<-predict(Alls.trans,newdata=Alls.preproc)
#dim(Alls.preproc2)

descrCor<-cor(Alls.preproc)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.85)
filteredDescr<-(Alls.preproc[,-highlyCorDescri])
dim(filteredDescr)

Asllls<-data.frame(filteredDescr, "bm"=Allstart[,bleachm])	
dim(Asllls)

set.seed(5)
boruta<-Boruta(bm~., data=Asllls, doTrace=F, maxRuns=500)
features<-attStats(boruta)

featuresConf<-rownames(subset(features, decision=="Confirmed"))

metrics<-c(featuresConf, "bm")		
dataNN<-Asllls[,match(metrics, colnames(Asllls))]
dim(dataNN)

sampsize<-list(20,25,30,35,40,45,50,55,60,65,70,75,80)
harvester<-list()
for(RR in 1:13)
{	#RR<-13
	resulter<-list()
	for(I in 1:100)
	{	
		set.seed(I)
		index = sample( seq_len ( nrow ( dataNN ) ), size = sampsize[[RR]] )
		
		trainNN<-dataNN[index,]
		testerNN<-dataNN[-index,]
		set.seed(I)
		indexer = sample( seq_len ( nrow ( testerNN ) ), size = 40 )
		testNN<-testerNN[indexer,]
		set.seed(I)
		model<-randomForest(bm ~., data=trainNN, mtry=1, ntree=4000)
		summary(model)	
		preds<-predict(model, testNN[,1:(dim(testNN)[2]-1)])
		
		modelEval<-cbind(testNN$bm, preds)
		colnames(modelEval)<-c('Actual','Predicted')
		modelEval<-as.data.frame(modelEval)
		
		R2<-cor(modelEval$Actual,modelEval$Predicted)^2
		Orep<-modelEval[order(modelEval$Predicted),]
		sigerBS<-data.frame(Orep[1:10,],"rep"=1)
		sigerBT<-data.frame(Orep[31:40,],"rep"=2)
		TvsS<-rbind(sigerBS, sigerBT)
		TvsSpval<-t.test(TvsS$Actual~TvsS$rep)[[3]]
		RMSE = (sum((modelEval$Actual - modelEval$Predicted)^2) / nrow(modelEval)) ^ 0.5
		BS<-mean(sigerBS$Actual)
		BT<-mean(sigerBT$Actual)
		
		resulter[[I]]<-data.frame("Tsize"=sampsize[[RR]], "rep"=I, "R2"=R2, "RMSE"=RMSE, "TvsS"=TvsSpval, "BSm"=BS, "BTm"=BT)
	}
	harvester[[RR]]<-do.call(rbind, resulter)
}
Harper<-do.call(rbind, harvester)

#write.csv(Harper, "AbsHarp1.csv")
Harper<-read.csv("AbsHarp1.csv", header=TRUE)

mRMSE<-aggregate(Harper$RMSE, by=list(Harper$Tsize), FUN=mean)
mRMSEsd<-aggregate(Harper$RMSE, by=list(Harper$Tsize), FUN=sd)
mr2<-aggregate(Harper$R2, by=list(Harper$Tsize), FUN=mean)
mr2sd<-aggregate(Harper$R2, by=list(Harper$Tsize), FUN=sd)
Harper1<-subset(Harper, TvsS<0.05)
mtest<-aggregate(Harper1$TvsS, by=list(Harper1$Tsize), FUN=NROW)
mRMSE
mr2
mtest



subset(Harper, Tsize=="80")
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

library(caret)
library(randomForest)
library(Boruta)

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

bleachm<-"Cha675"
Alls.cntn<-Allstart[,1:916]
Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]
dim(Alls.preproc)
Alls.trans<-preProcess(x=Alls.preproc,method=c("center","scale"))
Alls.preproc2<-predict(Alls.trans,newdata=Alls.preproc)
#dim(Alls.preproc2)

descrCor<-cor(Alls.preproc)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.85)
filteredDescr<-(Alls.preproc[,-highlyCorDescri])
dim(filteredDescr)

Asllls<-data.frame(filteredDescr, "bm"=Allstart[,bleachm])	
dim(Asllls)

set.seed(2)
boruta<-Boruta(bm~., data=Asllls, doTrace=F, maxRuns=500)
features<-attStats(boruta)
featuresConf<-rownames(subset(features, decision=="Confirmed"))

metrics<-c(featuresConf, "bm")		
dataNN<-Asllls[,match(metrics, colnames(Asllls))]
dim(dataNN)

I<-40
RR<-80
set.seed(I)
index = sample( seq_len ( nrow ( dataNN ) ), size = 80 )

trainNN<-dataNN[index,]
testerNN<-dataNN[-index,]
set.seed(I)
indexer = sample( seq_len ( nrow ( testerNN ) ), size = 40 )
testNN<-testerNN[indexer,]
set.seed(I)
model<-randomForest(bm ~., data=trainNN, mtry=1, ntree=4000)
summary(model)	
preds<-predict(model, testNN[,1:(dim(testNN)[2]-1)])

modelEval<-cbind(testNN$bm, preds)
colnames(modelEval)<-c('Actual','Predicted')
modelEval<-as.data.frame(modelEval)

R2<-cor(modelEval$Actual,modelEval$Predicted)^2
Orep<-modelEval[order(modelEval$Predicted),]
sigerBS<-data.frame(Orep[1:10,],"rep"=1)
sigerBT<-data.frame(Orep[31:40,],"rep"=2)
TvsS<-rbind(sigerBS, sigerBT)
TvsSpval<-t.test(TvsS$Actual~TvsS$rep)[[3]]
RMSE = (sum((modelEval$Actual - modelEval$Predicted)^2) / nrow(modelEval)) ^ 0.5
BS<-mean(sigerBS$Actual)
BT<-mean(sigerBT$Actual)

data.frame("Tsize"=RR, "rep"=I, "R2"=R2, "RMSE"=RMSE, "TvsS"=TvsSpval, "BSm"=BS, "BTm"=BT)

modelEval<-modelEval[order(modelEval$Predicted, decreasing=TRUE),]	


samp<-paste("MOTE-675Model", ".pdf", sep="")
pdf(file = samp, width = 10, height = 3, bg="transparent")

par(mfrow=c(1,3), oma = c(2, 2.5, 1.0, 2.5), mar = c(0.1, 3.1, 0.1, 5.1))

coll<-c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3)
plot(modelEval$Actual, modelEval$Predicted, ylim=c(-100,50), xlim=c(-100,50), axes=F, ylab="", xlab="", pch=19, col=c("black", "gray", "red")[coll])
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)
abline(lm(modelEval$Predicted~modelEval$Actual), col="red")

plot(mtest$Group.1, mtest$x, ylim=c(0,100), axes=F, pch=19, col="black", ylab="", xlim=c(20,80))
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)

rhose<-mr2sd$x/sqrt(100)
rhoT<-qt(p=0.05/2, df=99, lower.tail=F)*rhose

rmsese<-mRMSEsd$x/sqrt(100)
rmseT<-qt(p=0.05/2, df=99, lower.tail=F)*rmsese

plot(mr2$Group.1, mr2$x, ylim=c(0,0.5), axes=F, pch=19, col="black", ylab="", xlim=c(20,80))
segments(mr2$Group.1, mr2$x-(rhoT), mr2$Group.1, mr2$x+(rhoT), lwd=1, angle=90, length=0.05, code=3)
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)
par(new=TRUE)
plot(mRMSE$Group.1, mRMSE$x, ylim=c(0,50), axes=F, pch=19, col="gray", ylab="", xlim=c(20,80))
segments(mRMSE$Group.1, mRMSE$x-(rmseT), mRMSE$Group.1, mRMSE$x+(rmseT), lwd=1, angle=90, length=0.05, code=3)
axis(4, at=, col.axis="black", las=1, cex.axis=1.5)

dev.off()

summary(Harper)
meanQuant<-data.frame(t(aggregate(Harper[,7:8], by=list(Harper$Tsize), FUN=mean)[13,2:3]))
sdQuant<-aggregate(Harper[,7:8], by=list(Harper$Tsize), FUN=sd)[13,2:3]


samp<-paste("MOTE-675Bar", ".pdf", sep="")
pdf(file = samp, width = 4, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(2, 2.5, 1.0, 2.5), mar = c(0.1, 3.1, 0.1, 2.1))
sets<-barplot(meanQuant[,1], ylim=c(-100,50), col = c("red", "black"), axes=F)
arrows(sets, meanQuant[,1]-sdQuant[,1], sets, meanQuant[,1]+sdQuant[,1], lwd=1, angle=90, length=0.05, code=3)
axis(2, at=, col.axis="black", las=1)
box(col="black")

dev.off()

###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########
###########

library(randomForest)
library(Boruta)

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

bleachm<-"ChaQuant"
Alls.cntn<-Allstart[,1:916]
Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]
dim(Alls.preproc)
Alls.trans<-preProcess(x=Alls.preproc,method=c("center","scale"))
Alls.preproc2<-predict(Alls.trans,newdata=Alls.preproc)
#dim(Alls.preproc2)

descrCor<-cor(Alls.preproc)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.75)
filteredDescr<-(Alls.preproc[,-highlyCorDescri])
dim(filteredDescr)

Asllls<-data.frame(filteredDescr, "bm"=Allstart[,bleachm])	
dim(Asllls)

set.seed(2)  ###4
boruta<-Boruta(bm~., data=Asllls, doTrace=F, maxRuns=500)
features<-attStats(boruta)
featuresConf<-rownames(subset(features, decision=="Confirmed"))

metrics<-c(featuresConf, "bm")		
dataNN<-Asllls[,match(metrics, colnames(Asllls))]
dim(dataNN)

sampsize<-list(20,25,30,35,40,45,50,55,60,65,70,75,80)
harvester<-list()
for(RR in 1:13)
{	#RR<-13
	resulter<-list()
	for(I in 1:100)
	{	
		set.seed(I)
		index = sample( seq_len ( nrow ( dataNN ) ), size = sampsize[[RR]] )
		
		trainNN<-dataNN[index,]
		testerNN<-dataNN[-index,]
		set.seed(I)
		indexer = sample( seq_len ( nrow ( testerNN ) ), size = 40 )
		testNN<-testerNN[indexer,]
		set.seed(I)
		model<-randomForest(bm ~., data=trainNN, mtry=1, ntree=4000)
		summary(model)	
		preds<-predict(model, testNN[,1:(dim(testNN)[2]-1)])
		
		modelEval<-cbind(testNN$bm, preds)
		colnames(modelEval)<-c('Actual','Predicted')
		modelEval<-as.data.frame(modelEval)
		
		R2<-cor(modelEval$Actual,modelEval$Predicted)^2
		Orep<-modelEval[order(modelEval$Predicted),]
		sigerBS<-data.frame(Orep[1:10,],"rep"=1)
		sigerBT<-data.frame(Orep[31:40,],"rep"=2)
		TvsS<-rbind(sigerBS, sigerBT)
		TvsSpval<-t.test(TvsS$Actual~TvsS$rep)[[3]]
		RMSE = (sum((modelEval$Actual - modelEval$Predicted)^2) / nrow(modelEval)) ^ 0.5
		BS<-mean(sigerBS$Actual)
		BT<-mean(sigerBT$Actual)
		
		resulter[[I]]<-data.frame("Tsize"=sampsize[[RR]], "rep"=I, "R2"=R2, "RMSE"=RMSE, "TvsS"=TvsSpval, "BSm"=BS, "BTm"=BT)
	}
	harvester[[RR]]<-do.call(rbind, resulter)
}
Harper<-do.call(rbind, harvester)

#write.csv(Harper, "QuantHarp1.csv")
Harper<-read.csv("QuantHarp1.csv", header=TRUE)

mRMSE<-aggregate(Harper$RMSE, by=list(Harper$Tsize), FUN=mean)
mRMSEsd<-aggregate(Harper$RMSE, by=list(Harper$Tsize), FUN=sd)
mr2<-aggregate(Harper$R2, by=list(Harper$Tsize), FUN=mean)
mr2sd<-aggregate(Harper$R2, by=list(Harper$Tsize), FUN=sd)
Harper1<-subset(Harper, TvsS<0.05)
mtest<-aggregate(Harper1$TvsS, by=list(Harper1$Tsize), FUN=NROW)
mRMSE
mr2
mtest


subset(Harper, Tsize=="80")
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

bleachm<-"ChaQuant"
Alls.cntn<-Allstart[,1:916]
Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]
dim(Alls.preproc)
Alls.trans<-preProcess(x=Alls.preproc,method=c("center","scale"))
Alls.preproc2<-predict(Alls.trans,newdata=Alls.preproc)
#dim(Alls.preproc2)

descrCor<-cor(Alls.preproc)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.75)
filteredDescr<-(Alls.preproc[,-highlyCorDescri])
dim(filteredDescr)

Asllls<-data.frame(filteredDescr, "bm"=Allstart[,bleachm])	
dim(Asllls)

set.seed(2)
boruta<-Boruta(bm~., data=Asllls, doTrace=F, maxRuns=500)
features<-attStats(boruta)
featuresConf<-rownames(subset(features, decision=="Confirmed"))

metrics<-c(featuresConf, "bm")		
dataNN<-Asllls[,match(metrics, colnames(Asllls))]
dim(dataNN)

I<-55
RR<-80
set.seed(I)
index = sample( seq_len ( nrow ( dataNN ) ), size = 80 )

trainNN<-dataNN[index,]
testerNN<-dataNN[-index,]
set.seed(I)
indexer = sample( seq_len ( nrow ( testerNN ) ), size = 40 )
testNN<-testerNN[indexer,]
set.seed(I)
model<-randomForest(bm ~., data=trainNN, mtry=1, ntree=4000)
summary(model)	
preds<-predict(model, testNN[,1:(dim(testNN)[2]-1)])

modelEval<-cbind(testNN$bm, preds)
colnames(modelEval)<-c('Actual','Predicted')
modelEval<-as.data.frame(modelEval)

R2<-cor(modelEval$Actual,modelEval$Predicted)^2
Orep<-modelEval[order(modelEval$Predicted),]
sigerBS<-data.frame(Orep[1:10,],"rep"=1)
sigerBT<-data.frame(Orep[21:30,],"rep"=2)
TvsS<-rbind(sigerBS, sigerBT)
TvsSpval<-t.test(TvsS$Actual~TvsS$rep)[[3]]
RMSE = (sum((modelEval$Actual - modelEval$Predicted)^2) / nrow(modelEval)) ^ 0.5
BS<-mean(sigerBS$Actual)
BT<-mean(sigerBT$Actual)

data.frame("Tsize"=RR, "rep"=I, "R2"=R2, "RMSE"=RMSE, "TvsS"=TvsSpval, "BSm"=BS, "BTm"=BT)

modelEval<-modelEval[order(modelEval$Predicted, decreasing=TRUE),]			
			
samp<-paste("MOTE-QuantModel", ".pdf", sep="")
pdf(file = samp, width = 10, height = 3, bg="transparent")

par(mfrow=c(1,3), oma = c(2, 2.5, 1.0, 2.5), mar = c(0.1, 3.1, 0.1, 5.1))

coll<-c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3)
plot(modelEval$Actual, modelEval$Predicted, ylim=c(-100,100), xlim=c(-100,100), axes=F, ylab="", xlab="", pch=19, col=c("black", "gray", "red")[coll])
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)
abline(lm(modelEval$Predicted ~ modelEval$Actual), col="red")

plot(mtest$Group.1, mtest$x, ylim=c(0,100), axes=F, pch=19, col="black", ylab="", xlim=c(20,80))
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)


rhose<-mr2sd$x/sqrt(100)
rhoT<-qt(p=0.05/2, df=99, lower.tail=F)*rhose

rmsese<-mRMSEsd$x/sqrt(100)
rmseT<-qt(p=0.05/2, df=99, lower.tail=F)*rmsese

plot(mr2$Group.1, mr2$x, ylim=c(0,0.5), axes=F, pch=19, col="black", ylab="", xlim=c(20,80))
segments(mr2$Group.1, mr2$x-(rhoT), mr2$Group.1, mr2$x+(rhoT), lwd=1, angle=90, length=0.05, code=3)
box(col="black")
axis(1, at=, col.axis="black", las=1, cex.axis=1.5)
axis(2, at=, col.axis="black", las=1, cex.axis=1.5)
par(new=TRUE)
plot(mRMSE$Group.1, mRMSE$x, ylim=c(0,50), axes=F, pch=19, col="gray", ylab="", xlim=c(20,80))
segments(mRMSE$Group.1, mRMSE$x-(rmseT), mRMSE$Group.1, mRMSE$x+(rmseT), lwd=1, angle=90, length=0.05, code=3)
axis(4, at=, col.axis="black", las=1, cex.axis=1.5)

dev.off()

meanQuant<-data.frame(t(aggregate(Harper[,7:8], by=list(Harper$Tsize), FUN=mean)[13,2:3]))
sdQuant<-aggregate(Harper[,7:8], by=list(Harper$Tsize), FUN=sd)[13,2:3]

samp<-paste("MOTE-QuantBar", ".pdf", sep="")
pdf(file = samp, width = 4, height = 4, bg="transparent")

par(mfrow=c(1,1), oma = c(2, 2.5, 1.0, 2.5), mar = c(0.1, 3.1, 0.1, 2.1))
sets<-barplot(meanQuant[,1], ylim=c(-100,50), col = c("red", "black"), axes=F)
arrows(sets, meanQuant[,1]-sdQuant[,1], sets, meanQuant[,1]+sdQuant[,1], lwd=1, angle=90, length=0.05, code=3)
axis(2, at=, col.axis="black", las=1)
box(col="black")

dev.off()

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

library(caret)
library(randomForest)
library(Boruta)

lmp <- function (modelobject) {
    if (class(modelobject) != "lm") stop("Not an object of class 'lm' ")
    f <- summary(modelobject)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)
}

Allstart<-subset(Alls, Alls[,925]<400)
dim(Allstart)

bleachm<-"Cha675"
Alls.cntn<-Allstart[,1:916]
Alls.preproc<-subset(Alls.cntn, select= -c(qm.1.1, qm.1.0, qm.1.3, NPQ.1.0, NPQ.1.1, NPQ.1.3))
connectRmv<-subset(metafluor1, metric=="Connect")$uni
Alls.preproc<-Alls.preproc[,-which(colnames(Alls.preproc) %in% connectRmv)]
dim(Alls.preproc)
Alls.trans<-preProcess(x=Alls.preproc,method=c("center","scale"))
Alls.preproc2<-predict(Alls.trans,newdata=Alls.preproc)
#dim(Alls.preproc2)

descrCor<-cor(Alls.preproc2)
highlyCorDescri<-findCorrelation(descrCor, cutoff=0.99)
filteredDescr<-(Alls.preproc2[,-highlyCorDescri])
dim(filteredDescr)

Asllls<-data.frame(filteredDescr, "bm"=Allstart[,bleachm])
dim(Asllls)
lmcorQuant<-list()
rhocorQuant<-list()
for (R in 1:ncol(filteredDescr))
{	lmcorQuant[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=summary(lm(Asllls[,R] ~ Asllls[,"bm"]))[[9]], "pval"=lmp(lm(Asllls[,R] ~ Asllls[,"bm"])))
	rhocorQuant[[R]]<-data.frame("metrics"=colnames(Asllls[R]),"r2"=cor.test(Asllls[,R], Asllls[,"bm"], method=c("pearson"))[4][[1]], "pval"=cor.test(Asllls[,R], Asllls[,"bm"], method=c("pearson"))[3][[1]])
}
lmQuant<-data.frame(do.call(rbind, lmcorQuant))

resultQuant<-data.frame(metafluor1[match(lmQuant$metrics, metafluor1$uni),], lmQuant)
#resultQuant<-resultQuant[order(abs(resultQuant$r2), decreasing=TRUE),][1:100,]
resultQuant<-subset(resultQuant, pval<1)
dim(resultQuant)

metrics<-c(resultQuant$metrics)		
dataNN<-Asllls[,match(metrics, colnames(Asllls))]
dim(dataNN)

##########################
##########################
##########################
##########################
##########################
##########################
##########################

allerie2<-t(dataNN)
cb<-pvclust(allerie2, method.dist="manhattan", parallel=TRUE, method.hclust="ward.D", nboot=10000)
cbb<-cb

samp<-paste("DENDROeuclidean-wardD2", ".pdf", sep="")
pdf(file = samp, width = 6, height = 6, bg="transparent")
plot(cbb)
dev.off()

plot(cbb)

cbclusters<-cutree(cb, k=4)
cb <- color_branches(cb$hclust, k = 4, col = c("#f4a582", "#e31a1c", "#92c5de","#4393c3")) # add color to the lines
 
##########################
##########################
##########################
##########################
##########################
##########################

samp<-paste("MOTE-circularDendro-FULL", ".pdf", sep="")
pdf(file = samp, width = 2, height = 2, bg="transparent")
circlize_dendrogram(as.dendrogram(cb), labels_track_height = NA, dend_track_height = 0.9, cex=0.25, labels=F, lwd=2)
dev.off()

##########################
##########################
##########################
##########################
##########################
##########################

cbclust<-data.frame(cbclusters)
cb1<-subset(cbclust, cbclusters=="1") #grey
cb2<-subset(cbclust, cbclusters=="2") #red
cb3<-subset(cbclust, cbclusters=="3") #blue
cb4<-subset(cbclust, cbclusters=="4") #blue

phen1<-data.frame(Allstart[match(rownames(cb1),rownames(Allstart)),], "phen"="1")
dim(phen1)
phen2<-data.frame(Allstart[match(rownames(cb2),rownames(Allstart)),], "phen"="2")
dim(phen2)
phen3<-data.frame(Allstart[match(rownames(cb3),rownames(Allstart)),], "phen"="3")
dim(phen3)
phen4<-data.frame(Allstart[match(rownames(cb4),rownames(Allstart)),], "phen"="4")
dim(phen4)
analysisPhen<-rbind(phen1,phen2,phen3,phen4)

shapiro.test(analysisPhen$Cha675)
tm<-aov(analysisPhen$Cha675 ~ analysisPhen$phen)
summary(tm)
TukeyHSD(tm)

meanAbs<-aggregate(analysisPhen$Cha675, by=list(analysisPhen$phen), FUN=mean, na.rm=TRUE)
sdAbs<-aggregate(analysisPhen$Cha675, by=list(analysisPhen$phen), FUN=sd, na.rm=TRUE)
nrowAbs<-aggregate(analysisPhen$Cha675, by=list(analysisPhen$phen), FUN=NROW)
erAbs<-sdAbs$x/sqrt(nrowAbs$x)

samp<-paste("MOTEbar-abs", ".pdf", sep="")
pdf(file = samp, width = 5, height = 7, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 3.5, 0.1, 2.5), mar = c(1.1, 1.1, 1.1, 1.1))
sets<-barplot(meanAbs$x, ylim=c(-80,25), col = c("#e31a1c","#f4a582", "#4393c3", "#92c5de")[as.factor(meanAbs$Group.1)], axes=F)
arrows(sets, meanAbs$x-erAbs, sets, meanAbs$x+erAbs, lwd=1, angle=90, length=0.05, code=3)
axis(2, at=, col.axis="black", las=1, cex.axis=2.5)
box(col="black", lwd=2)
dev.off()

shapiro.test(sqrt(100+analysisPhen$ChaQuant))
tm<-aov(sqrt(100+analysisPhen$ChaQuant) ~ analysisPhen$phen)
summary(tm)
TukeyHSD(tm)

meanAbs<-aggregate(analysisPhen$ChaQuant, by=list(analysisPhen$phen), FUN=mean, na.rm=TRUE)
sdAbs<-aggregate(analysisPhen$ChaQuant, by=list(analysisPhen$phen), FUN=sd, na.rm=TRUE)
nrowAbs<-aggregate(analysisPhen$ChaQuant, by=list(analysisPhen$phen), FUN=NROW)
erAbs<-sdAbs$x/sqrt(nrowAbs$x)

samp<-paste("MOTEbar-absPred", ".pdf", sep="")
pdf(file = samp, width = 5, height = 7, bg="transparent")
par(mfrow=c(1,1), oma = c(2, 3.5, 0.1, 2.5), mar = c(1.1, 1.1, 1.1, 1.1))
sets<-barplot(meanAbs$x, ylim=c(-80,25), col = c("#e31a1c","#f4a582", "#4393c3", "#92c5de")[as.factor(meanAbs$Group.1)], axes=F)
arrows(sets, meanAbs$x-erAbs, sets, meanAbs$x+erAbs, lwd=1, angle=90, length=0.05, code=3)
#axis(2, at=, col.axis="black", las=1, cex.axis=2.5)
box(col="black", lwd=2)
dev.off()

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################

dim(FullPheno)
Ph1<-FullPheno[,1:1008][match(rownames(cb1), rownames(FullPheno)),]
Ph2<-FullPheno[,1:1008][match(rownames(cb2), rownames(FullPheno)),]
Ph3<-FullPheno[,1:1008][match(rownames(cb3), rownames(FullPheno)),]
Ph4<-FullPheno[,1:1008][match(rownames(cb4), rownames(FullPheno)),]

Phen1<-data.frame(metafluor1,"mean"=rowMeans(t(Ph1), na.rm=TRUE))
Phen2<-data.frame(metafluor1,"mean"=rowMeans(t(Ph2), na.rm=TRUE))
Phen3<-data.frame(metafluor1,"mean"=rowMeans(t(Ph3), na.rm=TRUE))
Phen4<-data.frame(metafluor1,"mean"=rowMeans(t(Ph4), na.rm=TRUE))

metris<-c("Quant", "qm", "qP", "NPQ", "Sigma", "ABQ", "Tau1", "Tau2")
sets<-c(0.4, 1, 1.2, 0.8, 12, 0.6, 6000, 4500000)
Phen<-Phen1

samp<-paste("Profiles-Phen1", ".pdf", sep="")
pdf(file = samp, width = 2, height = 7, bg="transparent")

par(mfrow=c(8,1), oma = c(2, 2, 1.5, 5), mar = c(0.5, 0.1, 0.5, 0.25))
for (X in 1:8)
{		
		col1<-subset(Phen, metric==metris[X] & col=="1")
		col2<-subset(Phen, metric==metris[X] & col=="2")
		col3<-subset(Phen, metric==metris[X] & col=="3")
		col4<-subset(Phen, metric==metris[X] & col=="4")
			
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)

		# if(X==1)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==2)
		# {	axis(4, at=seq(0,1.2,0.6), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==3)
		# {	axis(4, at=seq(0,0.8,0.4), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==4)
		# {	axis(4, at=seq(0,10,5), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==5)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==6)
		# {	axis(4, at=seq(0,0.6,0.3), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==7)
		# {	axis(4, at=seq(0,5000,2500), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==8)
		# {	axis(4, at=seq(0,4000000,2000000), col.axis="black", las=2, cex.axis=1.5)}
		if(X==8)
		{	axis(1, at=, col.axis="black", las=1, cex.axis=1.5)}
		box(col="black")
		grid(col="gray")

}
dev.off()


Phen<-Phen2

samp<-paste("Profiles-Phen2", ".pdf", sep="")
pdf(file = samp, width = 2, height = 7, bg="transparent")

par(mfrow=c(8,1), oma = c(2, 2, 1.5, 5), mar = c(0.5, 0.1, 0.5, 0.25))
for (X in 1:8)
{		
		col1<-subset(Phen, metric==metris[X] & col=="1")
		col2<-subset(Phen, metric==metris[X] & col=="2")
		col3<-subset(Phen, metric==metris[X] & col=="3")
		col4<-subset(Phen, metric==metris[X] & col=="4")
			
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)

		# if(X==1)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==2)
		# {	axis(4, at=seq(0,1.2,0.6), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==3)
		# {	axis(4, at=seq(0,0.8,0.4), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==4)
		# {	axis(4, at=seq(0,10,5), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==5)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==6)
		# {	axis(4, at=seq(0,0.6,0.3), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==7)
		# {	axis(4, at=seq(0,5000,2500), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==8)
		# {	axis(4, at=seq(0,4000000,2000000), col.axis="black", las=2, cex.axis=1.5)}
		if(X==8)
		{	axis(1, at=, col.axis="black", las=1, cex.axis=1.5)}
		box(col="black")
		grid(col="gray")

}
dev.off()


Phen<-Phen3

samp<-paste("Profiles-Phen3", ".pdf", sep="")
pdf(file = samp, width = 2, height = 7, bg="transparent")

par(mfrow=c(8,1), oma = c(2, 2, 1.5, 5), mar = c(0.5, 0.1, 0.5, 0.25))
for (X in 1:8)
{		
		col1<-subset(Phen, metric==metris[X] & col=="1")
		col2<-subset(Phen, metric==metris[X] & col=="2")
		col3<-subset(Phen, metric==metris[X] & col=="3")
		col4<-subset(Phen, metric==metris[X] & col=="4")
			
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)

		# if(X==1)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==2)
		# {	axis(4, at=seq(0,1.2,0.6), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==3)
		# {	axis(4, at=seq(0,0.8,0.4), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==4)
		# {	axis(4, at=seq(0,10,5), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==5)
		# {	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==6)
		# {	axis(4, at=seq(0,0.6,0.3), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==7)
		# {	axis(4, at=seq(0,5000,2500), col.axis="black", las=2, cex.axis=1.5)}
		# if(X==8)
		# {	axis(4, at=seq(0,4000000,2000000), col.axis="black", las=2, cex.axis=1.5)}
		if(X==8)
		{	axis(1, at=, col.axis="black", las=1, cex.axis=1.5)}
		box(col="black")
		grid(col="gray")

}
dev.off()


Phen<-Phen4

samp<-paste("Profiles-Phen4", ".pdf", sep="")
pdf(file = samp, width = 2, height = 7, bg="transparent")

par(mfrow=c(8,1), oma = c(2, 2, 1.5, 5), mar = c(0.5, 0.1, 0.5, 0.25))
for (X in 1:8)
{		
		col1<-subset(Phen, metric==metris[X] & col=="1")
		col2<-subset(Phen, metric==metris[X] & col=="2")
		col3<-subset(Phen, metric==metris[X] & col=="3")
		col4<-subset(Phen, metric==metris[X] & col=="4")
			
		plot(col1$stage, col1$mean, type="l", lwd=1.25, col="purple", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col2$stage, col2$mean, type="l", lwd=1.25, col="blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col3$stage, col3$mean, type="l", lwd=1.25, col="light blue", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)
		plot(col4$stage, col4$mean, type="l", lwd=1.25, col="cyan", ylim=c(0,sets[X]), axes=F)
		par(new=TRUE)

		if(X==1)
		{	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.25)}
		if(X==2)
		{	axis(4, at=seq(0,1,0.5), col.axis="black", las=2, cex.axis=1.25)}
		if(X==3)
		{	axis(4, at=seq(0,1.2,0.6), col.axis="black", las=2, cex.axis=1.25)}
		if(X==4)
		{	axis(4, at=seq(0,0.8,0.4), col.axis="black", las=2, cex.axis=1.25)}
		if(X==5)
		{	axis(4, at=seq(0,10,5), col.axis="black", las=2, cex.axis=1.25)}
		if(X==6)
		{	axis(4, at=seq(0,0.4,0.2), col.axis="black", las=2, cex.axis=1.25)}
		if(X==7)
		{	axis(4, at=seq(0,5000,2500), col.axis="black", las=2, cex.axis=1.25)}
		if(X==8)
		{	axis(4, at=seq(0,4000000,2000000), col.axis="black", las=2, cex.axis=1.25)}
		if(X==8)
		{	axis(1, at=, col.axis="black", las=1, cex.axis=1.5)}
		box(col="black")
		grid(col="gray")

}
dev.off()

##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
##########################
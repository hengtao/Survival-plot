library(survival)
library(survival)
library(GGally)
library(scales)
library(ggplot2)
setwd("C:/Users/Hengtao/Desktop")
(data<-read.table(file="lung_clinic_omic.txt",header=TRUE,na.strings = "NA",check.names = F))
survData<-data[,1:13]

data[,32]=factor(data[,32])
data$Subtype=data[,32]
data[,33]=factor(data[,33])
data$EGFR_NAV3=data[,33]
data[,11]=factor(data[,11])
data$ICDO_position=data[,11]
data[,14]=factor(data[,14])
data$EGFR=data[,14]

#单独绘制Subtype OS曲线图
fit<-coxph(Surv(data$time,as.numeric(data$Status)) ~ data$Subtype, data=data)
summary(fit)


fit <- survfit(Surv(as.numeric(survData$time),as.numeric(survData$Status))~data$Subtype,data=survData)
sub.text1<-"p=0.00124"
sub.text2<-"HR=0.1568"
sub.text3<-"n=95"
png(filename=paste(colnames(data)[32],"_os_cox.png",sep=""))
plot(fit,mark.time=TRUE,xlab="Time(Months)",ylab="Survival rate",main="Survival Curve",col=c("blue","red"),lty=c(2,3),lwd=2)
legend(6,0.10,legend=c(sub.text1),bty="n",col=1,cex=1.3)
legend("topright",c("Non-mutation","Mutation"),col=c("blue","red"),lty=c(2,3),lwd=2,cex=1.0)
legend(50,1.05,legend=sub.text3,bty="n",col=1,cex=1.5)
legend(40,.10,legend=sub.text2,bty="n",col=1,cex=1.3)
dev.off()


#单独绘制EGFR OS曲线图
fit<-coxph(Surv(data$time,as.numeric(data$Status)) ~ data$EGFR, data=data)
summary(fit)

fit <- survfit(Surv(as.numeric(survData$time),as.numeric(survData$Status))~data$EGFR,data=survData)
sub.text1<-"p=0.0500"
sub.text2<-"HR=2.1964"
sub.text3<-"n=95"
png(filename=paste(colnames(data)[14],"_os_cox.png",sep=""))
plot(fit,mark.time=TRUE,xlab="Time(Months)",ylab="Survival rate",main="Survival Curve",col=c("blue","red"),lty=c(2,3),lwd=2)
legend(6,0.10,legend=c(sub.text1),bty="n",col=1,cex=1.3)
legend("topright",c("Non-mutation","Mutation"),col=c("blue","red"),lty=c(2,3),lwd=2,cex=1.0)
legend(50,1.05,legend=sub.text3,bty="n",col=1,cex=1.5)
legend(40,.10,legend=sub.text2,bty="n",col=1,cex=1.3)
dev.off()

#单独绘制ICDO OS曲线图
fit<-coxph(Surv(data$time,as.numeric(data$Status)) ~ data$ICDO_position, data=data)
summary(fit)

fit <- survfit(Surv(as.numeric(survData$time),as.numeric(survData$Status))~data$ICDO_position,data=data)
sub.text1<-"p=0.0259"
sub.text2<-"HR=3.35517"
sub.text3<-"n=91"
sub.text4<-"p=0.9381"
sub.text5<-"HR=1.04291"
png(filename=paste(colnames(data)[11],"_os_cox.png",sep=""))
plot(fit,mark.time=TRUE,xlab="Time(Months)",ylab="Survival rate",col=c("blue","red","black"),lty=c(2,3),lwd=2)
legend(6,0.10,legend=c(sub.text1),bty="n",col=1,cex=1.3)
legend("topright",c("C34.1","C34.2,C34.3","C34.8,C34.9"),col=c("blue","red","black"),lty=c(2,3),lwd=2,cex=1.0)
legend(50,1.05,legend=sub.text3,bty="n",col=1,cex=1.5)
legend(40,.10,legend=sub.text2,bty="n",col=1,cex=1.3)
legend(6,0.18,legend=c(sub.text4),bty="n",col=1,cex=1.3)
legend(40,.18,legend=sub.text5,bty="n",col=1,cex=1.3)
dev.off()




#绘制该三个显著因素综合的OS曲线图
fit<-coxph(Surv(data$time,as.numeric(data$Status)) ~ data$EGFR_Subtype_ICDO, data=data)
summary(fit)

fit <- survfit(Surv(as.numeric(survData$time),as.numeric(survData$Status))~data$EGFR_Subtype_ICDO)
sub.text1<-"p=0.0053"
sub.text2<-"HR=4.0952"
sub.text3<-"n=91"
png(filename=paste(colnames(data)[34],"_os_Cox.png",sep=""))
plot(fit,mark.time=TRUE,xlab="Time(Months)",ylab="Survival rate",col=c("blue","red"),lty=c(2,3),lwd=2)
legend(6,0.10,legend=c(sub.text1),bty="n",col=1,cex=1.3)
legend("topright",c("EGFR=0,Subtype=1,ICDO!=1","Other"),col=c("blue","red"),lty=c(2,3),lwd=2,cex=1.0)
legend(50,1.05,legend=sub.text3,bty="n",col=1,cex=1.5)
legend(40,.10,legend=sub.text2,bty="n",col=1,cex=1.3)
dev.off()




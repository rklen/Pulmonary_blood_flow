linearInter<-function(x1,x,y){
  if(x1==x[1]){return(y[1])}
  for(j in 1:(length(x)-1)){
    if(x1>x[j] & x1<=x[j+1]){
      return(y[j]+(x1-x[j])/(x[j+1]-x[j])*(y[j+1]-y[j]))
    }
  }
}
fit_model<-function(param,input,j){
  k1<-abs(param[1])
  k2<-abs(param[2])
  Va<-1/(1+exp(-param[3]))
  Vb<-1/(1+exp(-param[4]))
  cT<-0
  for(i in 2:length(input)){
    if(i-j-1<1){c0<-0}else{
      if(i-j-1>length(input)){c0<-input[length(input)]}else{c0<-input[i-j-1]}
    }
    cT<-c(cT,k1*c0+(1-k2)*cT[i-1])
  }
  v<-(1-Va-Vb)*cT+Vb*input
  return(v)
}
errorfunction<-function(param,input,organCurve,j){
  return(sum((organCurve-fit_model(param,input,j))^2))
}
meanRelError<-function(modelCurve,organCurve){
  r<-0
  for(i in 2:length(organCurve)){
    if(modelCurve[i]!=organCurve[i]){
      r<-r+abs(modelCurve[i]-organCurve[i])/organCurve[i]
    }
  }
  r<-r/(length(organCurve)-1)
  return(r)
}
optimizationAlg<-function(initialValues,input,organCurve,j){
  tryCatch(
    {
      res<-nlm(errorfunction,initialValues,input=input,
               organCurve=organCurve,j=j,stepmax=1000,print.level=1)
      return(res$estimate)
    },
    error=function(e){
      return(initialValues)
    }
  )
}
colMedians<-function(d){
  medians<-c()
  for(i in 1:ncol(d)){
    medians[i]<-median(d[,i])
  }
  return(medians)
}
colSds<-function(d){
  sds<-c()
  for(i in 1:ncol(d)){
    sds[i]<-sd(d[,i])
  }
  return(sds)
}
wilcoxSymbol<-function(x,y){
  s<-''
  p<-wilcox.test(x,y,alternative='less',paired=TRUE)$p.value
  if(p<0.05){s<-'*'}
  if(p<0.01){s<-'**'}
  if(p<0.001){s<-'***'}
  return(s)
}

filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df<-read.table(filepath)
df<-df[c(2:9,11:28,30:34,36:61),]
sum(df$V3=='Female')

studyNumbers<-df$V1
t<-c(0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,80,90,100,120,140,
     160,190,220,250,280)
times<-c(0:280)

#1TCM with time delay, V_a, and V_b
v1<-c()
df1<-matrix(data=NA,nrow=length(studyNumbers),ncol=7)
for(i in 1:length(studyNumbers)){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/lun/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,7)
  input<-c()
  organCurve<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
    organCurve[l]<-linearInter(times[l],t,df[7,])
  }
  v1<-c(v1,max(input))
  #plot(times,input,type='l')
  plot(times,organCurve,type='l')
  initialValues<-c(1,1,-log(9),-log(9))
  errorsForj<-c()
  for(j in 0:15){
    param<-optimizationAlg(initialValues,input,organCurve,j)
    err<-errorfunction(param,input,organCurve,j)
    errorsForj<-c(errorsForj,err)
  }
  j<-c(0:15)[which.min(errorsForj)]
  #j<-0
  param<-optimizationAlg(initialValues,input,organCurve,j)
  points(times,fit_model(param,input,j),type='l')
  err<-errorfunction(param,input,organCurve,j)
  df1[i,1:length(param)]<-param
  df1[i,5]<-j
  df1[i,6]<-err/(length(times)-1)
  df1[i,7]<-meanRelError(fit_model(param,input,j),organCurve)
}
write.csv(df1,file=paste('lun_6.csv',sep=''),row.names=F)

df1<-read.csv('lun_6.csv')
filepath<-'C:/Users/Oona/Documents/Tpc/lok/Info for Oona.txt'
df<-read.table(filepath)
df<-df[2:61,]
df$V4<-c(90,94,79,85,73,82,81,108,63,92,64,117,99,59,89,100,103,
  101,78,86,87,77,114,81,90,78,65,94,121,61,95,64,60,92,
  101,69,62,64,132,98,130,67,123,80,NA,102,93,76,91,72,87,
  88,78,113,110,81,95,69,107,NA)
df$V5<-c(355,356,351,308,316,336,364,354,350,406,357,319,334,350,
  356,341,339,332,324,322,408,402,364,385,363,334,380,364,
  336,351,362,366,396,390,359,364,334,345,379,390,348,366,
  321,332,295,370,356,371,323,307,374,383,333,381,349,338,
  353,360,362,366)
df$V6<-c(168,180,172,169,180,172,170,173,166,174,158,185,164,166,
         159,177,172,170,179,176,174,178,184,160,164,167,163,164,
         176,152,183,161,158,192,180,164,158,161,167,171,182,167,
         190,171,NA,180,182,162,167,165,172,163,181,170,192,165,182,
         168,170,NA)
df$V7<-df[,4]/(df[,6]/100)^2
df<-df[c(1:8,10:27,29:33,35:60),]
df$V8<-abs(df1[,1])*60
df$V9<-abs(df1[,2])*60
df$V10<-abs(df1[,1]/df1[,2])
df$V11<-1/(1+exp(-df1[,3]))
df$V12<-1/(1+exp(-df1[,4]))
df$V13<-df1[,5]
df<-df[c(1:5,7:32,34:57),]
colMedians(df)
colMeans(df[,8:13])
colSds(df[,8:13])

plot(v1,df[,8])

womendf<-df[df[,3]=="Female",]
mendf<-df[df[,3]=="Male",]
wilcox.test(womendf[,8],mendf[,8],paired=FALSE)

colMedians(womendf)
colMedians(mendf)

youngwomendf<-womendf[womendf[,2]<65,]
oldwomendf<-womendf[womendf[,2]>64,]
colMedians(youngwomendf)
colMedians(oldwomendf)

youngmendf<-mendf[mendf[,2]<65,]
oldmendf<-mendf[mendf[,2]>64,]
colMedians(youngmendf)
colMedians(oldmendf)

youngdf<-df[df[,2]<65,]
olddf<-df[df[,2]>64,]
colMedians(youngdf)
colMedians(olddf)

wilcox.test(youngwomendf[,8],oldwomendf[,8],paired=FALSE)
wilcox.test(youngmendf[,8],oldmendf[,8],paired=FALSE)
wilcox.test(youngdf[,8],olddf[,8],paired=FALSE)
wilcox.test(youngwomendf[,8],youngmendf[,8],paired=FALSE)
wilcox.test(oldwomendf[,8],oldmendf[,8],paired=FALSE)
wilcox.test(womendf[,8],mendf[,8],paired=FALSE)

plot(as.numeric(mendf[,2]),mendf[,8],xlim=c(40,85),ylim=c(0,10),
     lwd=2,ylab='PBF (mL/min/cm^3)',xlab='Age (years)',cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[,2]),womendf[,8],pch=3,axes=FALSE,
       lwd=2)
fit<-lm(mendf[mendf[,8]<7,8]~as.numeric(mendf[mendf[,8]<7,2]))
abline(fit,col='blue',lw=2)
fit<-lm(womendf[womendf[,8]<7,8]~as.numeric(womendf[womendf[,8]<7,2]))
abline(fit,lw=2,lty=2)
legend('topright',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)
cor(as.numeric(mendf[mendf[,8]<7,2]),mendf[mendf[,8]<7,8])
cor(as.numeric(womendf[womendf[,8]<7,2]),womendf[womendf[,8]<7,8])

plot(as.numeric(mendf[!is.na(mendf[,4]),4]),
     mendf[!is.na(mendf[,4]),8],xlim=c(55,135),ylim=c(0,10),
     lwd=2,ylab='PBF (mL/min/cm^3)',xlab='Weight (kg)',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[!is.na(womendf[,4]),4]),
       womendf[!is.na(womendf[,4]),8],pch=3,axes=FALSE,
       lwd=2)
fit<-lm(mendf[!is.na(mendf[,4])&mendf[,8]<7,8]~
          as.numeric(mendf[!is.na(mendf[,4])&mendf[,8]<7,4]))
abline(fit,col='blue',lw=2)
fit<-lm(womendf[!is.na(womendf[,4])&womendf[,8]<7,8]~
          as.numeric(womendf[!is.na(womendf[,4])&womendf[,8]<7,4]))
abline(fit,lw=2,lty=2)
legend('topright',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)
cor(as.numeric(mendf[!is.na(mendf[,4])&mendf[,8]<7,4]),
    mendf[!is.na(mendf[,4])&mendf[,8]<7,8])
cor(as.numeric(womendf[!is.na(womendf[,4])&womendf[,8]<7,4]),
    womendf[!is.na(womendf[,4])&womendf[,8]<7,8])

plot(as.numeric(mendf[!is.na(mendf[,7]),7]),
     mendf[!is.na(mendf[,7]),8],xlim=c(20,50),ylim=c(0,10),
     lwd=2,ylab='PBF (mL/min/cm^3)',xlab='BMI',
     cex.axis=1.3,cex.lab=1.3,col='blue')
points(as.numeric(womendf[!is.na(womendf[,7]),7]),
       womendf[!is.na(womendf[,7]),8],pch=3,axes=FALSE,
       lwd=2)
fit<-lm(mendf[!is.na(mendf[,7])&mendf[,8]<7,8]~
          as.numeric(mendf[!is.na(mendf[,7])&mendf[,8]<7,7]))
abline(fit,col='blue',lw=2)
fit<-lm(womendf[!is.na(womendf[,7])&womendf[,8]<7,8]~
          as.numeric(womendf[!is.na(womendf[,7])&womendf[,8]<7,7]))
abline(fit,lw=2,lty=2)
legend('topright',legend=c('Women','Men'),lty=c(2,1),pch=c(3,1),
       col=c('black','blue'),lwd=c(2,2),cex=1.4)
cor(as.numeric(mendf[!is.na(mendf[,7])&mendf[,8]<7,7]),
    mendf[!is.na(mendf[,7])&mendf[,8]<7,8])
cor(as.numeric(womendf[!is.na(womendf[,7])&womendf[,8]<7,7]),
    womendf[!is.na(womendf[,7])&womendf[,8]<7,8])

cor(as.numeric(df[df[,8]<7,2]),df[df[,8]<7,8])
cor(as.numeric(df[!is.na(df[,4])&df[,8]<7,4]),
    df[!is.na(df[,4])&df[,8]<7,8])
cor(as.numeric(df[!is.na(df[,7])&df[,8]<7,7]),
    df[!is.na(df[,7])&df[,8]<7,8])

mean(as.numeric(df[,2]))
sd(as.numeric(df[,2]))
min(as.numeric(df[,2]))
max(as.numeric(df[,2]))

min(as.numeric(df[,5]))
max(as.numeric(df[,5]))

filepath<-paste('C:/Users/Oona/Documents/Tpc/lun/array1_',studyNumbers[1],'.csv',sep='')
df<-read.csv(filepath)
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,7)
input<-c()
organCurve<-c()
organCurve1<-c()
organCurve2<-c()
organCurve3<-c()
organCurve4<-c()
organCurve5<-c()
for(l in 1:length(times)){
  #input[l]<-linearInter(times[l],t,df[1,])
  organCurve1[l]<-linearInter(times[l],t,df[2,])
  organCurve2[l]<-linearInter(times[l],t,df[3,])
  organCurve3[l]<-linearInter(times[l],t,df[4,])
  organCurve4[l]<-linearInter(times[l],t,df[5,])
  organCurve5[l]<-linearInter(times[l],t,df[6,])
  organCurve[l]<-linearInter(times[l],t,df[7,])
}
plot(times,organCurve5/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
     col='blue',lty=2,cex.axis=1.3,cex.lab=1.3,
     ylab='Acitivity concentration (MBq/mL)',xlab='Time (s)')
points(times,organCurve1/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
       col='gray')
points(times,organCurve2/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
       col='gray',lty=2)
points(times,organCurve3/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
       col='steelblue1',lty=2)
points(times,organCurve4/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
       col='blue')
points(times,organCurve/1000,type='l',lwd=2,ylim=c(0,45),xlim=c(0,280),
       col='black')
legend('topright',legend=c('LLL','RLL','RML','LUL','RUL','All lobes'),
       lty=c(1,2,2,1,2,1),col=c('gray','gray','steelblue1','blue',
                                'blue','black'),lwd=c(2,2,2,2,2,2),
       cex=1.3)

filepath<-paste('C:/Users/Oona/Documents/Tpc/lun/array1_',studyNumbers[1],'.csv',sep='')
df<-read.csv(filepath)
df<-as.matrix(df)
colnames(df)<-NULL
rownames(df)<-NULL
df[,1]<-rep(0,7)
input<-c()
for(l in 1:length(times)){
  input[l]<-linearInter(times[l],t,df[1,])
}
plot(times,input/1000,type='l',lwd=2,ylim=c(0,230),xlim=c(0,280),
     cex.axis=1.3,cex.lab=1.3,
     ylab='Acitivity concentration (MBq/mL)',xlab='Time (s)')
for(i in 2:10){
  filepath<-paste('C:/Users/Oona/Documents/Tpc/lun/array1_',studyNumbers[i],'.csv',sep='')
  df<-read.csv(filepath)
  df<-as.matrix(df)
  colnames(df)<-NULL
  rownames(df)<-NULL
  df[,1]<-rep(0,7)
  input<-c()
  for(l in 1:length(times)){
    input[l]<-linearInter(times[l],t,df[1,])
  }
  points(times,input/1000,type='l',lwd=2,ylim=c(0,230),xlim=c(0,280),
       axis=FALSE)
}

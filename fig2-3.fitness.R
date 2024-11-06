library(plyr)
library(dplyr)
library(parallel)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
###filter raw data
logphase <- function(x){
  tmp <- combn(nrow(x),2) %>% t() %>% as.data.frame() %>% dplyr::filter(V2-V1>1)
  tmpRP <- lapply(1:nrow(tmp), function(y){
    a <- x[tmp[y,1]:tmp[y,2],]
    t <- summary(lm(OD~time,data = a))
    data.frame(stringsAsFactors = F,DTslope=log(2)/data.frame(t$coefficients)[2,1],p=data.frame(t$coefficients)[2,4],Rsquare=t$r.squared) %>%
      cbind(a)
  }) %>% rbind.fill() %>% dplyr::filter(p<0.05) 
  mybest <- (tmpRP[,2:3] %>% unique() %>% arrange(desc(Rsquare)))[1:ceiling(nrow(tmp)*0.4),] %>% dplyr::filter(p==min(p))
  tmpRP %>% dplyr::filter(p==mybest$p & Rsquare==mybest$Rsquare)
  
}
allFiles <- system("ls ~/project/HGT/exp/fitness.epoach/GmR_*",intern=T)
alldata <- mclapply(mc.cores=10,1:length(allFiles),function(x) {
  a <- as.data.frame(read.csv(allFiles[x], skip=66,nrows=100,sep = "\t",fileEncoding="latin1",fill=T))[,c(-1,-2)]
  type <- as.numeric(strsplit(allFiles[x],"_")[[1]][2])
  a$time <- seq(from=1,to=1000,by=10)
  blank <- 0.089
  colrun <- colnames(a)[which((substr(colnames(a),1,3)== "X0."))]
  lapply(1:length(colrun), function(xx){
    data <- data.frame(OD=a[,which(colnames(a)==colrun[xx])]-blank)
    data$time <- a$time
    data$name <- colrun[xx]
    data$dp <- as.numeric(strsplit(strsplit(colrun[xx],c("_"))[[1]][1],"X")[[1]][2])
    data$seq <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][2])
    data$biorep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][3])
    data$tecrep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][4])
    data$batch <- strsplit(colrun[xx],c("_"))[[1]][5]
    if(max(data$OD)>1){
      
      data %>% dplyr::filter(OD > 0.2 & OD < 0.6) -> mydf
      mydf$OD <- log(mydf$OD)
      #log phase   
      LP <- logphase(mydf)
      
      DTvector <- log(2)/((sort(diff(LP$OD)/diff(LP$time),decreasing = T)))
      DTvector <- DTvector[which(DTvector != Inf)]
      res <- mydf[1,3:8] %>% cbind(data.frame(stringsAsFactors = F,type,DTmin1=DTvector[1],DTmin2=DTvector[2],
                                              DTmean=mean(DTvector),DTmedian=median(DTvector),DTslope=unique(LP$DTslope)))
      
    }else {
      res <- data[1,3:8] %>% cbind(data.frame(stringsAsFactors = F,type,DTmin1=99,DTmin2=99,
                                              DTmean=99,DTmedian=99,
                                              DTslope=99))
    }
    
    res
  }) %>% rbind.fill()
}) %>% rbind.fill() %>% dplyr::filter(DTmean!=99)
alldata %>% dplyr::filter(dp!=0) -> alldata
save(alldata,file = "~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.Rdata")

### filter with cv
options(dplyr.summarise.inform = FALSE)
source("~/Rfunction/style.print.R")
load("~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.Rdata")
mydata0 <- alldata %>% dplyr::filter(dp!=0) %>%
  group_by(batch,type,dp,seq,biorep) %>%
  dplyr::mutate(cvfinalDTmin1=sd(DTmin1)/mean(DTmin1)) %>%
  as.data.frame()
a <- mydata0[,c(2:7,8,13)]
names(a)[7:8] <- c("m","cv")
cvcutoff1 <- 0.15
df0 <- a %>% dplyr::filter(cv<cvcutoff1)
df1 <- a %>% dplyr::filter(cv>=cvcutoff1) %>% dplyr::mutate(zz=paste(batch,dp,seq,biorep,type)) %>% group_by(zz) %>% dplyr::filter(length(zz)>=3)
df2 <- mclapply(mc.cores = 10,1:length(unique(df1$zz)), function(x){
  mya <- df1 %>% dplyr::filter(zz==unique(df1$zz)[x]) %>% as.data.frame()
  tmp1 <- lapply(1:nrow(mya), function(y){
    b <- mya[y,]
    c <- mya[-y,]
    data.frame(stringsAsFactors = F,y,dis=sum(abs(b$m-c$m)))
  }) %>% rbind.fill() %>% arrange(dis)
  newa <- mya[c(tmp1$y[1:2]),] %>% dplyr::mutate(cv=sd(m)/mean(m))
  if(unique(newa$cv)<cvcutoff1){
    res <- newa[,1:8]
  } else{
    res <- data.frame(stringsAsFactors = F,cv=999) %>% cbind(df0[1,1:7])
  }
  res
  
}) %>% rbind.fill() %>% dplyr::filter(cv!=999)

mydf <- rbind(df0,df2)
WT <- (mydf %>% dplyr::filter(type == 0))$m %>% mean()
mydf$fitness0 <- WT/mydf$m
save(mydf,file = "~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.filter.Rdata")

###plot
##cut dp with three groups
load("~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.filter.Rdata")
load("~/project/HGT/exp/newdp.Rdata")
DP %>% ggplot(aes(x=oldDP,y=newDP))+geom_point()+geom_vline(xintercept = c(1/3,2/3))+geom_hline(yintercept = c(1/3,0.6))

mydf$newdp <- DP$newDP[match(mydf$dp,DP$oldDP)]
names(mydf)[c(1,10)] <- c("olddp","dp")
highlow.P <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x]) %>% group_by(dp,seq,biorep) %>% dplyr::summarize(fitness0=mean(fitness0))
  rank1 <- a %>% dplyr::filter(dp < 1/3)
  rank2 <- a %>% dplyr::filter(dp > 1/3 & dp < 2/3)
  rank3 <- a %>% dplyr::filter(dp > 2/3)
  data.frame(stringsAsFactors = F,
             m=c(mean(rank1$fitness0)/mean(rank3$fitness0),mean(rank2$fitness0)/mean(rank3$fitness0),mean(rank3$fitness0)/mean(rank3$fitness0)),
             se=c(sd(rank1$fitness0/mean(rank3$fitness0))/(length(rank1$fitness0)^0.5),sd(rank2$fitness0/mean(rank3$fitness0))/(length(rank2$fitness0)^0.5),sd(rank3$fitness0/mean(rank3$fitness0))/(length(rank3$fitness0)^0.5)),
             typevalue=c("(0,1/3)","(1/3,2/3)","(2/3,1)"),
             type=c("rank1","rank2","rank3"),
             typeGmR=unique(mydf$type)[x])
  
}) %>% rbind.fill() %>% arrange(typeGmR)
highlow.P$typeGmR <- as.character(highlow.P$typeGmR)
highlow.P$typeGmR <- factor(highlow.P$typeGmR,levels=c("0","10","20","30","40","60","70","80","100","120","140","150","160","180"))
highlow.P$typevalue <- factor(highlow.P$typevalue,levels=c("(0,1/3)","(1/3,2/3)","(2/3,1)"))
source("~/Rfunction/style.print.R")
highlow.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(typeGmR,m,fill=typevalue)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=m+se,ymin=m-se,color=typevalue),position = position_dodge(width = 0.8),width = 0.3)+
  #scale_fill_manual(values = c(rep("#E6E6E6",2),rep("#B4B4B4",8),rep("#4D4D4D",3)))+
  #scale_y_continuous(limits = c(0,1.5),breaks = c(0,0.5,1,1.5))+
  labs(x="Gentamicin concentration (ug/ml)",y="Relative fitness (relative to the large mismatch group)")+
  style.print()

### p values
options(scipen = 200)
highlow.Pvalue <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x]) %>% group_by(dp,seq,biorep) %>% dplyr::summarize(fitness0=mean(fitness0))
  rank1 <- a %>% dplyr::filter(dp < 1/3)
  rank2 <- a %>% dplyr::filter(dp > 1/3 & dp < 2/3)
  rank3 <- a %>% dplyr::filter(dp > 2/3)
  p12 <- wilcox.test(rank1$fitness0,rank2$fitness0)$p.value
  p13 <- wilcox.test(rank1$fitness0,rank3$fitness0)$p.value
  p23 <- wilcox.test(rank2$fitness0,rank3$fitness0)$p.value
  data.frame(stringsAsFactors = F,typeGmR=unique(mydf$type)[x],p12,p13,p23)
}) %>% rbind.fill() %>% arrange(typeGmR)


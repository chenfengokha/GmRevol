library(plyr)
library(dplyr)
library(parallel)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
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


allFiles <- system("ls /home/chenfeng/project/HGT/exp/01.amp/02.fitness/amp*txt",intern=T)

alldata <- mclapply(mc.cores=10,1:length(allFiles),function(x) {
  a <- as.data.frame(read.csv(allFiles[x], skip=83,nrows=100,sep = "\t",fileEncoding="latin1",fill=T))[,c(-1,-2)]
  type <- as.numeric(strsplit(allFiles[x],"-")[[1]][3])
  day <- as.numeric(substr(strsplit(allFiles[x],"-")[[1]][4],4,4))
  a$time <- seq(from=1,to=1000,by=10)
  blank <- 0.08772333
  #a[,which(substr(colnames(a),1,2)=="LB")] %>% colMeans() %>% mean()
  colrun <- colnames(a)[which((substr(colnames(a),1,3) %in% c("X0.","WT_")))]
  lapply(1:length(colrun), function(xx){
    data <- data.frame(OD=a[,which(colnames(a)==colrun[xx])]-blank)
    data$time <- a$time
    data$name <- colrun[xx]
    data$dp <- as.numeric(strsplit(strsplit(colrun[xx],c("_"))[[1]][1],"X")[[1]][2])
    data$biorep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][2])
    data$tecrep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][3])
   
    if(max(data$OD)>0.8){
      
      data %>% dplyr::filter(OD > 0.2 & OD < 0.6) -> mydf
      mydf$OD <- log(mydf$OD)
        
      LP <- logphase(mydf)
      
      DTvector <- log(2)/((sort(diff(LP$OD)/diff(LP$time),decreasing = T)))
      DTvector <- DTvector[which(DTvector != Inf)]
      res <- mydf[1,3:6] %>% cbind(data.frame(stringsAsFactors = F,day,type,DTmin1=DTvector[1],DTmin2=DTvector[2],
                                              DTmean=mean(DTvector),DTmedian=median(DTvector),DTslope=unique(LP$DTslope)))
    }else {
      res <- data[1,3:6] %>% cbind(data.frame(stringsAsFactors = F,day,type,DTmin1=99,DTmin2=99,
                                              DTmean=99,DTmedian=99,
                                              DTslope=99))
    }
    res
  }) %>% rbind.fill()
  
}) %>% rbind.fill() %>% dplyr::filter(DTmean!=99)

tmpcontrol <- alldata %>% dplyr::filter(is.na(dp)) %>% group_by(day,type,biorep) %>% dplyr::summarize(Dmin1=mean(DTmin1),Dmin2=mean(DTmin2),Dmean=mean(DTmean),Dmedian=mean(DTmedian),Dslope=mean(DTslope))
experimentdata <- lapply(1:length(unique(tmpcontrol$day)), function(x){
  experimentdf <- alldata %>% dplyr::filter(!is.na(dp) & day==unique(tmpcontrol$day)[x])

  controldf <- tmpcontrol %>% dplyr::filter(day==unique(tmpcontrol$day)[x])
  experimentdf$Dmin1.1 <- controldf$Dmin1[match(experimentdf$type,controldf$type)]
  experimentdf$Dmin2.1 <- controldf$Dmin2[match(experimentdf$type,controldf$type)]
  experimentdf$Dmean.1 <- controldf$Dmean[match(experimentdf$type,controldf$type)]
  experimentdf$Dmedian.1 <- controldf$Dmedian[match(experimentdf$type,controldf$type)]
  experimentdf$Dslope.1 <- controldf$Dslope[match(experimentdf$type,controldf$type)]

  experimentdf %>%
    dplyr::filter(!is.na(Dmin1.1)) %>%
    group_by(day,type,dp,biorep,tecrep) %>%
    dplyr::summarize(Dmin1=(DTmin1/Dmin1.1),
                     Dmin2 =(DTmin2/Dmin2.1),
                     Dmean =(DTmean/Dmean.1),
                     Dmedian=(DTmedian/Dmedian.1),
                     Dslope =(DTslope/Dslope.1))

}) %>% rbind.fill()

WT <- (experimentdata %>% dplyr::filter(type == 100))$Dmin1 %>% mean()
experimentdata$fitness0 <- WT/experimentdata$Dmin1
load("/home/chenfeng/project/HGT/exp/01.amp/02.fitness/amp.DP.Rdata")
experimentdata$DP0.1 <- testdp$DP0.1[match(experimentdata$dp,round(testdp$oldDP,3))]
names(experimentdata)[c(3,12)] <- c("olddp","dp")

experimentdata %>% group_by(type) %>%
  dplyr::summarize(rho=cor.test(dp,fitness0)$estimate,
                   P=cor.test(dp,fitness0)$p.value)
experimentdata %>% ggplot(aes(dp,fitness0))+geom_point()+geom_smooth(method = "lm",se=F)+facet_grid(.~type)

#######filter data

mydata0 <- experimentdata %>% dplyr::filter(dp!=0) %>%
  group_by(type,biorep,dp) %>%
  dplyr::mutate(cvfinalDTmin1=sd(Dmin1)/mean(Dmin1)) %>%
  as.data.frame()

a <- mydata0[,c(1,2,12,4:6,13)]

names(a)[6:7] <- c("m","cv")

cvcutoff1 <- 0.20
df0 <- a %>% dplyr::filter(cv<cvcutoff1)
df1 <- a %>% dplyr::filter(cv>=cvcutoff1) %>% dplyr::mutate(zz=paste(dp,biorep,type)) %>% group_by(zz) %>% dplyr::filter(length(zz)>=3)

df2 <- mclapply(mc.cores = 10,1:length(unique(df1$zz)), function(x){
  mya <- df1 %>% dplyr::filter(zz==unique(df1$zz)[x]) %>% as.data.frame()
  tmp1 <- lapply(1:nrow(mya), function(y){
    b <- mya[y,]
    c <- mya[-y,]
    data.frame(stringsAsFactors = F,y,dis=sum(abs(b$m-c$m)))
  }) %>% rbind.fill() %>% arrange(dis)
  (mya[c(tmp1$y[1:2]),] %>% dplyr::mutate(cv=sd(m)/mean(m)))[,-8]
  
}) %>% rbind.fill() %>% dplyr::filter(cv<cvcutoff1)

mydf <- rbind(df0,df2)

WT <- (mydf %>% dplyr::filter(type == 100))$m %>% mean()
mydf$fitness0 <- WT/mydf$m
mydf <- mydf %>% group_by(dp,biorep,type) %>% dplyr::summarize(m=mean(m),fitness0=mean(fitness0))
mydf %>% group_by(type) %>%
  dplyr::summarize(rho=cor.test(dp,fitness0)$estimate,
                   P=cor.test(dp,fitness0)$p.value)
save(mydf,file = "/home/chenfeng/project/HGT/exp/01.amp/02.fitness/alldata.all.RPbest.filter.AmpR.Rdata")


##########################################
load("/home/chenfeng/project/HGT/exp/01.amp/02.fitness/alldata.all.RPbest.filter.AmpR.Rdata")
source("~/Rfunction/style.print.R")
tmpaa <- data.frame(stringsAsFactors = F,dp=unique(mydf$dp) %>% sort(),rank=1:length(unique(mydf$dp)))
###########Fig. S5a
tmpaa %>% ggplot(aes(x=rank,y=dp))+geom_point()+scale_x_continuous(limits = c(1,47),breaks = c(16,32))+scale_y_continuous(breaks = c(0.2,0.4,0.6))+style.print()
###############test of Fig.S5b S 
mydf %>% group_by(type) %>%
  dplyr::summarize(rho=cor.test(dp,fitness0,method="s")$estimate,
                   P=cor.test(dp,fitness0,method="s")$p.value)


#################################

source("~/Rfunction/style.print.R")
highlow.P <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x])
  
  rank1 <- a %>% dplyr::filter(dp <= 0.375)
  rank2 <- a %>% dplyr::filter(dp > 0.375 & dp <= 0.485)
  rank3 <- a %>% dplyr::filter(dp > 0.485)
  
  data.frame(stringsAsFactors = F,
             len=c(length(rank1$dp %>% unique()),length(rank2$dp %>% unique()),length(rank3$dp %>% unique())),
             m=c(mean(rank1$fitness0)/mean(rank3$fitness0),mean(rank2$fitness0)/mean(rank3$fitness0),mean(rank3$fitness0)/mean(rank3$fitness0)),
             se=c(sd(rank1$fitness0/mean(rank3$fitness0))/(length(rank1$fitness0)^0.5),sd(rank2$fitness0/mean(rank3$fitness0))/(length(rank2$fitness0)^0.5),sd(rank3$fitness0/mean(rank3$fitness0))/(length(rank3$fitness0)^0.5)),
             type=c("rank1","rank2","rank3"),
             typeGmR=rep(unique(mydf$type)[x],3))
  
}) %>% rbind.fill() %>% arrange(typeGmR)
highlow.P$typeGmR <- as.character(highlow.P$typeGmR)
highlow.P$typeGmR <- factor(highlow.P$typeGmR,levels=c("100","2000","4000"))
highlow.P$type <- factor(highlow.P$type,levels=c("rank1","rank2","rank3","rank4","rank5"))
source("~/Rfunction/style.print.R")
########Fig. S5b
highlow.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(typeGmR,m,fill=type)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=m+se,ymin=m-se),position = position_dodge(width = 0.8),width = 0.3)+
  #scale_fill_manual(values = c(rep("#E6E6E6",2),rep("#B4B4B4",8),rep("#4D4D4D",3)))+
  scale_y_continuous(limits = c(0,2),breaks = c(0,1,2))+
  labs(x="Ampicillin concentration (ug/ml)",y="Relative fitness (relative to the large mismatch group)")+
  style.print()

### test of Fig. S5b
options(scipen = 200)
highlow.Pvalue <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x])
  rank1 <- a %>% dplyr::filter(dp <= 0.375)
  rank2 <- a %>% dplyr::filter(dp > 0.375 & dp <= 0.485)
  rank3 <- a %>% dplyr::filter(dp > 0.485)
  p12 <- wilcox.test(rank1$fitness0,rank2$fitness0)$p.value
  p13 <- wilcox.test(rank1$fitness0,rank3$fitness0)$p.value
  
  p23 <- wilcox.test(rank2$fitness0,rank3$fitness0)$p.value
  
  data.frame(stringsAsFactors = F,typeGmR=unique(mydf$type)[x],p12,p13,p23)
}) %>% rbind.fill() %>% arrange(typeGmR)


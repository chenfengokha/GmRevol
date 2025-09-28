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


allFiles <- system("ls /home/chenfeng/project/HGT/exp/Gmr_tRNA/newdata/*.txt",intern=T)

alldata <- mclapply(mc.cores=3,1:length(allFiles),function(x) {
  a <- as.data.frame(read.csv(allFiles[x], skip=84,nrows=145,sep = "\t",fileEncoding="latin1",fill=T))[,c(-1,-2)]
  a$time <- seq(from=1,to=1450,by=10)
  blank <- 0.089
  colrun <- colnames(a)[which((substr(colnames(a),1,3)== "X0."))]
  tt <- lapply(1:length(colrun), function(xx){
    data <- data.frame(OD=a[,which(colnames(a)==colrun[xx])]-blank)
    data$time <- a$time
    data$name <- colrun[xx]
    data$dp <- as.numeric(strsplit(strsplit(colrun[xx],c("_"))[[1]][1],"X")[[1]][2])
    data$seq <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][2])
    data$biorep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][3])
    data$tecrep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][4])
    data$GmRtype <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][5])
    data$trnatype <- strsplit(colrun[xx],c("_"))[[1]][6]
    if(max(data$OD)>0.6){
      
      data %>% dplyr::filter(OD > 0.2 & OD < 0.6) -> mydf
      mydf$OD <- log(mydf$OD)
      #log phase   
      LP <- logphase(mydf)
      
      DTvector <- log(2)/((sort(diff(LP$OD)/diff(LP$time),decreasing = T)))
      DTvector <- DTvector[which(DTvector != Inf)]
      res <- mydf[1,3:9] %>% cbind(data.frame(stringsAsFactors = F,DTmin1=DTvector[1],DTmin2=DTvector[2],
                                              DTmean=mean(DTvector),DTmedian=median(DTvector),DTslope=unique(LP$DTslope)))
    }else {
      res <- data[1,3:9] %>% cbind(data.frame(stringsAsFactors = F,DTmin1=99,DTmin2=99,
                                              DTmean=99,DTmedian=99,
                                              DTslope=99))
    }
    if(unique(res$trnatype)=="n"){
      res$type2 <- "Orignal"
    } else{
      res$type2 <- "Novel"
    }
    res
    
  }) %>% rbind.fill()
}) %>% rbind.fill() %>% dplyr::filter(DTmean!=99)


##filter data
mydf <- alldata
mydf$fitness <- mean(mydf$DTmin1)/mydf$DTmin1

mydf$zz <- paste(mydf$dp,mydf$seq,mydf$biorep,mydf$GmRtype,mydf$type2)

fdata <- mclapply(mc.cores = 10,1:length(unique(mydf$zz)), function(x){
  a <- (mydf %>% dplyr::filter(zz == unique(mydf$zz)[x]) %>% dplyr::mutate(delta=abs(median(DTmin1)-DTmin1)) %>% arrange(delta))
  if(nrow(a)>9){
    res <- a[1:9,]
  } else {
    res <- a
  }
  res
}) %>% rbind.fill()

fdata$type3 <- ""
fdata$type3[which(fdata$dp==0.778 & fdata$type2=="Novel")] <- "tRNA (argW) supply increased, D=0.632 (small)"
fdata$type3[which(fdata$dp==0.778 & fdata$type2=="Orignal")] <- "Wildtype, D=0.778 (large)"
fdata$type3[which(fdata$dp==0.328 & fdata$type2=="Novel")] <- "tRNA (argW) supply increased, D=0.338 (small)"
fdata$type3[which(fdata$dp==0.328 & fdata$type2=="Orignal")] <- "Wildtype, D=0.328 (minimum)" 

save(fdata,file = "/home/chenfeng/project/HGT/exp/Gmr_tRNA/new/merge.fitness.GmR.tRNA.Rdata")
##############plot
load("/home/chenfeng/project/HGT/exp/Gmr_tRNA/new/merge.fitness.GmR.tRNA.Rdata")
source("~/Rfunction/style.print.R")
fdata$GmRtype <- as.character(fdata$GmRtype)
fdata$GmRtype <- factor(fdata$GmRtype,levels = c("10","120","180"))
fdata$type3 <- factor(fdata$type3,levels = c("Wildtype, D=0.328 (minimum)","tRNA (argW) supply increased, D=0.338 (small)","tRNA (argW) supply increased, D=0.632 (small)","Wildtype, D=0.778 (large)"))
fdata$GmRtypedp <- paste(fdata$GmRtype,fdata$dp)
finaldata <- lapply(1:length(unique(fdata$GmRtypedp)), function(x){
  mya <- fdata %>% dplyr::filter(GmRtypedp == unique(fdata$GmRtypedp)[x])
  nn <- mya$fitness[which(mya$trnatype == "argw")] %>% mean()
  mya %>% group_by(GmRtype,type2,type3,dp,seq) %>%
    dplyr::summarize(m=mean(fitness/nn),se=sd(fitness/nn)/(length(fitness)^0.5))
}) %>% rbind.fill()

##Fig. S4c
finaldata %>%
  ggplot(aes(x=GmRtype,y=m,fill=type3))+ 
  geom_bar(stat = "identity",position = position_dodge(preserve = 'single'))+
  geom_errorbar(aes(ymax=m+se,ymin=m-se),width=0.2,position = position_dodge(0.9))+
  scale_fill_manual(values=c("#F8766D","#00BA38","#00BA38","#619CFF"))+
  facet_grid(dp~.,scales = "free")+
  scale_y_continuous(limits = c(0,2),breaks = c(0,1,2))+
  labs(x="Gentamicin concentration (μg/ml)",y="Fitness")+
  style.print()

##test of Fig. S4c
##0.778
a1 <- fdata %>% dplyr::filter(dp==0.778 & seq==1 & type2=="Novel" & GmRtype ==10)
a2 <- fdata %>% dplyr::filter(dp==0.778 & seq==1 & type2=="Orignal" & GmRtype ==10)
wilcox.test(a1$fitness,a2$fitness,alternative = "less")

b1 <- fdata %>% dplyr::filter(dp==0.778 & seq==1 & type2=="Novel" & GmRtype ==120)
b2 <- fdata %>% dplyr::filter(dp==0.778 & seq==1 & type2=="Orignal" & GmRtype ==120)
wilcox.test(b1$fitness,b2$fitness,alternative = "greater")

c1 <- fdata %>% dplyr::filter(dp==0.778 & seq==1 & type2=="Novel" & GmRtype ==180)
wilcox.test(c1$fitness,mu=0, alternative = "greater")

##0.328
a3 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Novel" & GmRtype ==10)
a4 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Orignal" & GmRtype ==10)
wilcox.test(a3$fitness,a4$fitness,alternative = "greater")

b3 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Novel" & GmRtype ==120)
b4 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Orignal" & GmRtype ==120)
wilcox.test(b3$fitness,b4$fitness,alternative = "greater")

c3 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Novel" & GmRtype ==180)
c4 <- fdata %>% dplyr::filter(dp==0.328 & seq==1 & type2=="Orignal" & GmRtype ==180)
wilcox.test(c3$fitness,c4$fitness, alternative = "less")



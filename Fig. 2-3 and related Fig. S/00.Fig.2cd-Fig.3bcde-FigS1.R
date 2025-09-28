options(dplyr.summarise.inform = FALSE)
library(magrittr);
library(parallel);
library(plyr);
library(dplyr);
library(gridExtra);
library(ggplot2);
library(flowCore);
library(lattice);
library(flowViz);
library(ggpmisc);

allFiles <- system("ls /mnt/data4/disk/chenfeng/HGT/exp/facs/data/GM80/*fcs",intern=T)
alldata <- mclapply(mc.cores=20,1:length(allFiles),function(x){
  a <- as.data.frame(exprs(read.FCS(allFiles[x],transformation="linearize", emptyValue = F)))
  b <- a[,c(1,2,4,9)] %>% dplyr::filter(Time>0)
  names(b) <- c("time","sizeA","YFPA","mCherryA")
  tmpdd <- strsplit(strsplit(strsplit(allFiles[x], "/")[[1]][11],".fcs")[[1]][1],"_")[[1]]
  #tmpdd <- strsplit(strsplit(strsplit(allFiles[x], "/")[[1]][11],".fcs")[[1]][1],"_")[[1]][-1]
  b$Dp <- as.numeric(tmpdd[5])
  b$strain <- as.numeric(tmpdd[6])
  b$biorep <- as.numeric(tmpdd[7])
  b$techrep <- as.numeric(tmpdd[8])
  b$GmCon <- as.numeric(tmpdd[9])
  b$batch <- tmpdd[10]
  b$type <- strsplit(strsplit(allFiles[x], "Group_")[[1]][2],".fcs")[[1]][1]
  b$type1 <- paste(tmpdd[5],tmpdd[6],tmpdd[7],tmpdd[9],sep = "_")
  
  mydf <- b %>% as.data.frame(stringsAsFactors = F) %>%
    dplyr::filter(Dp > 0 &  sizeA < 100000 & sizeA > 10000) %>%
    dplyr::mutate(YFP=YFPA/sizeA,mCherry=mCherryA/sizeA) %>%
    arrange(mCherry)
  mydf1 <- mydf[ceiling(nrow(mydf)*0.025):ceiling(nrow(mydf)*0.975),] %>% arrange(YFP)
  mydf1$type2 <- paste(mydf1$GmCon,mydf1$batch)
  
  c <- mydf1[ceiling(nrow(mydf1)*0.025):ceiling(nrow(mydf1)*0.975),] %>% 
    dplyr::filter(YFP > 10*0.01370595 &  mCherry > 10*0.0001204383) %>%
    arrange(mCherry)
  if(nrow(c)/nrow(b)>0.2){
    res <- c
  } else{res <- cbind(c[1,-15],data.frame(stringsAsFactors = F,type2="999"))}
  
  res
  
}) %>% rbind.fill() %>% dplyr::filter(type2!="999")


########1.0  0,10,20,30,40,60,80,100,120,140,160,180 
###merge data
load("~/project/HGT/exp/facs/data/merge_0.Rdata")
df0 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP),
                   Size=mean(sizeA),
                   mexp=mean(sizeA*mCherryA),
                   yexp=mean(sizeA*YFPA)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP),
                   smSize=mean(Size),
                   smexp=mean(mexp),
                   syexp=mean(yexp)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)

load("~/project/HGT/exp/facs/data/merge_10.Rdata")
df10 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP),
                   Size=mean(sizeA),
                   mexp=mean(mCherryA),
                   yexp=mean(YFPA)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP),
                   smSize=mean(Size),
                   smexp=mean(mexp),
                   syexp=mean(yexp)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)

load("~/project/HGT/exp/facs/data/merge_20.Rdata")
df20 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_30.Rdata")
df30 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_40.Rdata")
df40 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_60.Rdata")
df60 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_80.Rdata")
df80 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_100.Rdata")
df100 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_120.Rdata")
df120 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_140.Rdata")
df140 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_160.Rdata")
df160 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)
load("~/project/HGT/exp/facs/data/merge_180.Rdata")
df180 <- alldata %>%
  group_by(batch,Dp,strain,biorep,techrep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA0=mean(mCherry),
                   smeanYFP=mean(YFP)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::mutate(cvmCherry=sd(SmeanmCherrydensityA0)/mean(SmeanmCherrydensityA0)) %>%
  group_by(batch,Dp,strain,GmCon) %>%
  dplyr::filter(cvmCherry==min(cvmCherry)) %>%
  group_by(batch,Dp,strain,biorep,GmCon) %>%
  dplyr::summarize(SmeanmCherrydensityA1=mean(SmeanmCherrydensityA0),
                   smeanYFP1=mean(smeanYFP)) %>%
  as.data.frame() %>% 
  dplyr::mutate(maxMcherry=max(SmeanmCherrydensityA1),
                maxYFP=max(smeanYFP1)) %>%
  dplyr::mutate(payoff=SmeanmCherrydensityA1/maxMcherry,
                cost=maxYFP/smeanYFP1)

alldata <- rbind(df0,df10,df20,df30,df40,df60,df80,df100,df120,df140,df160,df180)
save(alldata,file="/mnt/data4/disk/chenfeng/HGT/exp/facs/data/draw.expression3.Rdata")

#############################################################
###3)plot
##xij estimated using CUB of top 10% highly expressed genes
load("/mnt/data4/disk/chenfeng/HGT/exp/facs/data/draw.expression3.Rdata")
names(alldata)[c(2)] <- c("dp")

##xij estimated using tRNA expression level
load("~/project/HGT/exp/fig2-3.related/gmr.Dtrnaexp.Rdata")
alldata$d.trnaexp <- newD.tRNAexp$DP.trnaexp[match(alldata$dp,newD.tRNAexp$oldDP)]
names(alldata)[c(2,12)] <- c("olddp","dp")

group.P <- lapply(1:length(unique(alldata$GmCon)), function(x){
  a <- alldata %>% dplyr::filter(GmCon==unique(alldata$GmCon)[x])
  rank1 <- a %>% dplyr::filter(dp < 0.33)
  rank2 <- a %>% dplyr::filter(dp > 0.33 & dp < 0.66)
  rank3 <- a %>% dplyr::filter(dp > 0.66)
  
  data.frame(stringsAsFactors = F,
             mLpayoffratio=rep(mean(rank1$payoff)/mean(rank3$payoff),3),
             mLcostratio=rep(mean(rank1$cost)/mean(rank3$cost),3),
             mPayoff=c(mean(rank1$payoff)/mean(rank3$payoff),mean(rank2$payoff)/mean(rank3$payoff),mean(rank3$payoff)/mean(rank3$payoff)),
             sePayoff=c(sd(rank1$payoff/mean(rank3$payoff))/(length(rank1$payoff)^0.5),sd(rank2$payoff/mean(rank3$payoff))/(length(rank2$payoff)^0.5),sd(rank3$payoff/mean(rank3$payoff))/(length(rank3$payoff)^0.5)),
             mCost=c(mean(rank1$cost)/mean(rank3$cost),mean(rank2$cost)/mean(rank3$cost),mean(rank3$cost)/mean(rank3$cost)),
             seCost=c(sd(rank1$cost/mean(rank3$cost))/(length(rank1$cost)^0.5),sd(rank2$cost/mean(rank3$cost))/(length(rank2$cost)^0.5),sd(rank3$cost/mean(rank3$cost))/(length(rank3$cost)^0.5)),
             typevalue=c("(0,1/3)","(1/3,2/3)","(2/3,1)"),
             type=c("rank1","rank2","rank3"),
             typeGmR=unique(alldata$GmCon)[x])
  
}) %>% rbind.fill() %>% arrange(typeGmR)
group.P$typeGmR <- as.character(group.P$typeGmR)
group.P$typeGmR <- factor(group.P$typeGmR,levels=c("0","10","20","30","40","60","70","80","100","120","140","150","160","180"))
group.P$typevalue <- factor(group.P$typevalue,levels=c("(0,1/3)","(1/3,2/3)","(2/3,1)"))
source("~/Rfunction/style.print.R")
##Fig.2c and Fig.3b
group.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(typeGmR,mPayoff,fill=typevalue)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=mPayoff+sePayoff,ymin=mPayoff-sePayoff,color=typevalue),position = position_dodge(width = 0.8),width = 0.3)+
  scale_y_continuous(limits = c(0,5),breaks = c(0,1,2,3,4))+
  labs(x="Gentamicin concentration (ug/ml)",y="Functional payoff (relative to the large mismatch group)")+
  style.print()
##Fig.2d and Fig.3d
group.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(typeGmR,mCost,fill=typevalue)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=mCost+seCost,ymin=mCost-seCost,color=typevalue),position = position_dodge(width = 0.8),width = 0.3)+
  scale_y_continuous(limits = c(0,2.5),breaks = c(0,1,2))+
  labs(x="Gentamicin concentration (ug/ml)",y="Translational cost (relative to the large mismatch group)")+
  style.print()
##wilcox.test of Fig.2c and Fig.3b
options(scipen = 200)
group.Pvalue <- lapply(1:length(unique(alldata$GmCon)), function(x){
  a <- alldata %>% dplyr::filter(GmCon==unique(alldata$GmCon)[x])
  rank1 <- a %>% dplyr::filter(dp < 0.33)
  rank2 <- a %>% dplyr::filter(dp > 0.33 & dp < 0.66)
  rank3 <- a %>% dplyr::filter(dp > 0.66)
  p12 <- wilcox.test(rank1$payoff,rank2$payoff)$p.value
  p13 <- wilcox.test(rank1$payoff,rank3$payoff)$p.value
  p23 <- wilcox.test(rank2$payoff,rank3$payoff)$p.value
  data.frame(stringsAsFactors = F,typeGmR=unique(alldata$GmCon)[x],p12,p13,p23)
}) %>% rbind.fill() %>% arrange(typeGmR)
##wilcox.test of Fig.2d and Fig.3d
group.Pvalue <- lapply(1:length(unique(alldata$GmCon)), function(x){
  a <- alldata %>% dplyr::filter(GmCon==unique(alldata$GmCon)[x])
  rank1 <- a %>% dplyr::filter(dp < 0.33)
  rank2 <- a %>% dplyr::filter(dp > 0.33 & dp < 0.66)
  rank3 <- a %>% dplyr::filter(dp > 0.66)
  p12 <- wilcox.test(rank1$cost,rank2$cost)$p.value
  p13 <- wilcox.test(rank1$cost,rank3$cost)$p.value
  p23 <- wilcox.test(rank2$cost,rank3$cost)$p.value
  data.frame(stringsAsFactors = F,typeGmR=unique(alldata$GmCon)[x],p12,p13,p23)
}) %>% rbind.fill() %>% arrange(typeGmR)


##Fig.3c,e
fig.3ce <- group.P %>% dplyr::filter(type=="rank1") 
fig.3ce$typeGmR <- as.numeric(as.vector(fig.3ce$typeGmR))
fig.3ce %>% ggplot(aes(x=typeGmR,y=mPayoff))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_y_continuous(limits = c(1,3.6),breaks = c(1,2,3))+
  labs(x="Gentamicin concentration (μg/ml)",y="Functional payoff of the minimum mismatch group\nrelative to the large mismatch group")+
  style.print()
fig.3ce %>% ggplot(aes(x=typeGmR,y=mCost))+
  geom_point()+
  geom_smooth(method = "lm",se=F)+
  scale_y_continuous(limits = c(1,2),breaks = c(1,1.5,2))+
  labs(x="Gentamicin concentration (μg/ml)",y="Translational cost of the minimum mismatch group\nrelative to the large mismatch group")+
  style.print()

cor.test(fig.3ce$typeGmR,fig.3ce$mPayoff,method = "s")
cor.test(fig.3ce$typeGmR,fig.3ce$mCost,method = "s")


###figs1ab
##tRNA exp
load("/mnt/data4/disk/chenfeng/HGT/exp/facs/data/draw.expression3.Rdata")

metmpdd <- "s"
payoff <- alldata %>% group_by(GmCon) %>%
  dplyr::summarize(rho=as.vector(cor.test(Dp,payoff,method = metmpdd )$estimate),
                   p=as.vector(cor.test(Dp,payoff,method = metmpdd )$p.value),
                   type="Functional payoff")
cost <- alldata %>% group_by(GmCon) %>%
  dplyr::summarize(rho=as.vector(cor.test(Dp,cost,method = metmpdd )$estimate),
                   p=as.vector(cor.test(Dp,cost,method = metmpdd )$p.value),
                   type="Translational cost")


###figs1c
load("~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.filter.Rdata")

fit <- lapply(1:length(unique(mydf$type)),function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x])
  
  b <- a %>%
    group_by(type,dp,seq,biorep) %>%
    dplyr::summarize(nn=mean(fitness0))
  if(as.numeric(cor.test(b$dp,b$nn,method = metmpdd)$estimate[1])>0){
    res2 <- data.frame(rho=as.vector(cor.test(b$dp,b$nn,method = metmpdd )$estimate[1]),
                       p=as.vector(cor.test(b$dp,b$nn,method = metmpdd )$p.value),
                       stringsAsFactors = F,GmCon=unique(mydf$type)[x],type="Fitness")
    
  } else{
    res2 <- data.frame(rho=as.vector(cor.test(b$dp,b$nn,method = metmpdd )$estimate[1]),
                       p=as.vector(cor.test(b$dp,b$nn,method = metmpdd )$p.value),
                       stringsAsFactors = F,GmCon=unique(mydf$type)[x],type="Fitness")
  }
  res2
  
}) %>% rbind.fill() %>% arrange(GmCon)
  
figs1 <- rbind(payoff,cost,fit)
figs1$type <- factor(figs1$type,levels = c("Fitness","Translational cost","Functional payoff"))
figs1$GmCon <- factor(figs1$GmCon,levels = unique(figs1$GmCon))
figs1$text <- NA
figs1$text[which(figs1$p<0.1)] <- "*"
figs1$text[which(figs1$p<0.05)] <- "**"
figs1$text[which(figs1$p<0.01)] <- "***"

figs1$rho[which(figs1$rho < 0)] <- figs1$rho[which(figs1$rho < 0)]*0.2/max(figs1$rho[which(figs1$rho < 0)] %>% abs())
figs1$rho[which(figs1$rho > 0)] <- figs1$rho[which(figs1$rho > 0)]*0.2/max(figs1$rho[which(figs1$rho > 0)] %>% abs())

source("~/Rfunction/style.print.R")
ggplot(figs1, aes(GmCon,type,fill=rho))+
  geom_tile(color="grey",alpha=2)+
  #geom_raster()+
  #scale_fill_manual(values = c("white","#0000FF","#FF0000"))+
  scale_fill_gradient2(low="navy",high = "firebrick3",mid="white")+
  geom_text(aes(label=text),col="black")+
  style.print()+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1, colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        plot.margin=unit(c(0.4,0.4,0.4,0.4),units=,"cm"))

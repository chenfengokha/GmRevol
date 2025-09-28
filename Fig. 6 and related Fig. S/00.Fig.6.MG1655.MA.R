library(plyr)
library(dplyr)
library(ggplot2)
library(parallel)
#library(Biostrings)
library(stringr)
library(tidyr)

#1)calculate DP for all gene with one mutation (simulation)
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
dpofref <- mclapply(mc.cores=6,1:nrow(nameseq),function(x){
  myname <- nameseq$bname[x]
  myseq <- as.character(nameseq$myseq[x])
  mutseq <- mclapply(mc.cores=5,4:(nchar(myseq)-3),function(y){
    seqtmp <- myseq
    str_sub(seqtmp,y,y) <- "A"
    aseq <- seqtmp
    str_sub(seqtmp,y,y) <- "T"
    bseq <- seqtmp
    str_sub(seqtmp,y,y) <- "G"
    cseq <- seqtmp
    str_sub(seqtmp,y,y) <- "C"
    dseq <- seqtmp
    data.frame(stringsAsFactors = F,seqall=c(aseq,bseq,cseq,dseq))
  }) %>% rbind.fill() %>% unique()
  mutseq$type <- "mut"
  mutseq$type[which(mutseq$seqall==myseq)] <- "wt"
  alldp <- mclapply(mc.cores=5,1:nrow(mutseq),function(z){
    datatype <- mutseq[z,2]
    dataseq <- mutseq[z,1]
    a <- substring(dataseq,seq(1,(nchar(dataseq)-2),by=3),seq(3,nchar(dataseq),by=3)) %>% table() %>% as.data.frame()
    names(a)[1] <- c("codon")
    a$codon <- as.vector(a$codon)
    a$amino <- codon$aa[match(a$codon,codon$codon)]
    yima <- sum((a %>% dplyr::filter(codon %in% c("TGA","TAG","TAA")))$Freq)
    if(yima==1){
      b <- a %>% dplyr::filter(amino != "NA" & amino != "W") %>%
        group_by(amino) %>% dplyr::mutate(nn=sum(Freq)) %>% dplyr::mutate(yij=Freq/nn)
      xijdata <- xijofD
      xijdata$yij <- b$yij[match(xijdata$codon,b$codon)]
      xijdata$yij[which(is.na(xijdata$yij))] <- 0
      res <- xijdata %>% group_by(amino) %>% 
        dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditai=(sum((yij-xijoftai)^2))^0.5,Ditaiexp=(sum((yij-xijoftRNAexp)^2))^0.5) %>%
        dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptai=prod(Ditai)^(1/length(Ditai)),Dptaiexp=prod(Ditaiexp)^(1/length(Ditaiexp))) %>%
        cbind(data.frame(stringsAsFactors = F,datatype))
    } else {
      res <- data.frame(stringsAsFactors = F,Dpexp=99,Dptai=99,datatype)
    }
    res
  }) %>% rbind.fill() %>% dplyr::filter(Dpexp!=99)
  if(length(unique(alldp$datatype))>1){
    finres <- alldp %>% dplyr::filter(datatype=="mut") %>% 
      cbind(data.frame(stringsAsFactors = F,WTdpexp=alldp$Dpexp[which(alldp$datatype=="wt")],wtdptai=alldp$Dptai[which(alldp$datatype=="wt")],wtdptaiexp=alldp$Dptaiexp[which(alldp$datatype=="wt")]))
    names(finres)[1:3] <- c("mutdpexp","mutdptai","mutdptaiexp")
  } else {
    finres <- data.frame(mutdpexp=99,mutdptai=99,datatype=99,WTdpexp=99,wtdptai=99,mutdptaiexp=99,wtdptaiexp=99,stringsAsFactors = F)
  }
  finres$bname <- myname 
  finres
}) %>% rbind.fill() %>% dplyr::filter(mutdptai!=99)

save(dpofref,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.ofrefgene.randommutation.mg1655.1.Rdata")
# # #2)calculate DP for MA genes
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
MAmutant1 <- read.csv("~/project/HGT/exp/model/MG1655MA.csv",stringsAsFactors = F) %>% dplyr::filter(type %in% c("WT","MMR") & Amino.acid.mutation=="-" & bname %in% nameseq$bname)
dpofmut <- mclapply(mc.cores=10,1:nrow(MAmutant1),function(x){
  atmp <- MAmutant1[x,]
  myfrom <- substr(atmp$Codon.mutation,1,3)
  myto <- substr(atmp$Codon.mutation,6,8)
  btmp <- nameseq %>% dplyr::filter(toupper(bname)==toupper(atmp$bname)) %>% dplyr::mutate(len=nchar(myseq)) %>% dplyr::filter(len==max(len))
  dataseq <- btmp$myseq
  codonofgene <- substring(dataseq,seq(1,(nchar(dataseq)-2),by=3),seq(3,nchar(dataseq),by=3)) %>% table() %>% as.data.frame()
  names(codonofgene)[1] <- "codon"
  codonofgene$codon <- as.vector(codonofgene$codon)
  codonofgene %>% rbind(data.frame(stringsAsFactors = F,codon=c(myfrom,myto),Freq=c(-1,1))) -> codonofgene
  codonofgene$aa <- codon$aa[match(codonofgene$codon,codon$codon)]
  yijMA <- codonofgene %>% dplyr::filter(!(aa %in% c(NA,"M","W"))) %>%
    group_by(aa,codon) %>% dplyr::summarize(newfre=sum(Freq)) %>%
    group_by(aa) %>% dplyr::mutate(aafre=sum(newfre)) %>%
    as.data.frame() %>% dplyr::mutate(yij=newfre/aafre)
  xijdata <- xijofD
  xijdata$yij <- yijMA$yij[match(xijdata$codon,yijMA$codon)]
  xijdata$yij[which(is.na(xijdata$yij))] <- 0
  xijdata %>% group_by(amino) %>%
    dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditai=(sum((yij-xijoftai)^2))^0.5,Ditaiexp=(sum((yij-xijoftRNAexp)^2))^0.5) %>%
    dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptai=prod(Ditai)^(1/length(Ditai)),Dptaiexp=prod(Ditaiexp)^(1/length(Ditaiexp))) %>%
    cbind(data.frame(stringsAsFactors = F,gene=atmp$Gene,bname=atmp$bname,pos=atmp$Position))

}) %>% rbind.fill()
save(dpofmut,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpofmutgene.MAline.MG1655.1.Rdata")

# #2)
##1)exp0.1 of D
load("/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.ofrefgene.randommutation.mg1655.1.Rdata")
dpofref <- dpofref[,c(1,5,8)]
names(dpofref)[1:2] <- c("mutdp","wtdp")
load("/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpofmutgene.MAline.MG1655.1.Rdata")
dpofmut <- dpofmut[,c(1,5)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$bname)]
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2

express <- express %>%
  arrange(desc(average)) %>% dplyr::filter(log.TPM %in% dpofref$bname & log.TPM %in% dpofmut$bname)

i=0.1
##111)top 10% highly expressed genes
experimentG <- express$log.TPM[1:floor(i*nrow(express))]

amutinf <- dpofmut %>% dplyr::filter(bname %in% experimentG) %>% group_by(bname) %>% dplyr::summarize(n=length(bname)) %>% group_by(n) %>% dplyr::summarize(nn=length(n))
amutinf1 <- dpofmut %>% dplyr::filter(bname %in% experimentG) %>% group_by(bname) %>% dplyr::summarize(n=length(bname))
mutif <- dpofmut %>% group_by(bname) %>% dplyr::summarize(nn=length(mutdp))

controlG <- lapply(1:length(amutinf$n), function(j){

  ainf <- amutinf[j,]
  ainf1 <- amutinf1 %>% dplyr::filter(n==ainf$n)
  binf <- mutif %>% dplyr::filter(nn==ainf$n & !(mutif$bname %in% experimentG))
  if(nrow(binf)==0){
    res <- data.frame(bname=ainf1$bname,nn=NA,Dp=NA)
  } else {
    binf$Dp <- dpofmut$wtdp[match(binf$bname,dpofmut$bname)]
    res <- (binf %>% arrange(Dp))[1:ainf$nn,]
  }
  res %>% cbind(type="Control")
}) %>% rbind.fill()

tmpgene <- data.frame(stringsAsFactors = F,type="Experiment",gene=setdiff(experimentG,controlG$bname)) %>%
  rbind(data.frame(stringsAsFactors = F,type="Control",gene=setdiff(controlG$bname,experimentG)))

k=1
typeoftwo <- unique(tmpgene$type)[k]
genetmp <- tmpgene$gene[which(tmpgene$type==typeoftwo)]

amut <- dpofmut %>% dplyr::filter(bname %in% genetmp)
bref <- dpofref %>% dplyr::filter(bname %in% genetmp)

atmp <- mclapply(mc.cores=30,1:length(genetmp), function(y){
  amutdata <- amut %>% dplyr::filter(bname==genetmp[y])
  brefall <- bref %>% dplyr::filter(bname==genetmp[y])

  subdata <- mclapply(mc.cores = 2,1:10000,function(x){
    set.seed(x)
    brefall[sample(nrow(brefall),nrow(amutdata),replace = T),] %>% cbind(data.frame(stringsAsFactors = F,x,type="con"))

  }) %>% rbind.fill()

  subdata %>% rbind(cbind(amutdata,data.frame(stringsAsFactors = F,x=10001,type="mut")))

}) %>% rbind.fill()

muttmp <- atmp %>% dplyr::filter(type=="mut")
contmp <- atmp %>% dplyr::filter(type=="con")

#1)
MArise <- (muttmp %>% dplyr::filter(mutdp>wtdp)  %>% nrow())/nrow(muttmp)
dataresrise <- contmp %>%
  group_by(x) %>%
  dplyr::mutate(nn=length(x)) %>% as.data.frame() %>%
  group_by(x) %>%
  dplyr::summarize(CONrise=length(which(mutdp>wtdp))/unique(nn))

source("~/Rfunction/style.print.R")
dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
  geom_segment(aes(x=MArise,y=200,xend=MArise,yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
  #scale_x_continuous(breaks = c(0.35,0.45,0.55))+
  labs(x="Percentage of D-increasing mutations \namong top 10% highly expressed genes",y="Frequency")+
  #labs(x="Percentage of D−increasing mutations\namong genes with 10% lowest D values",y="Frequency")+
  style.print()
data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > MArise)/10000,typeoftwo,i)
sd(dataresrise$CONrise)
mean(dataresrise$CONrise)

#2)
MArise <- (muttmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0))$deltadp
dataresrise <- contmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0) %>%
  group_by(x) %>%
  dplyr::summarize(CONrise=mean(deltadp))

source("~/Rfunction/style.print.R")
dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
  geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
  #scale_x_continuous(breaks = c(5,15,25))+
  labs(x="Average effect sizes of D-increasing mutations",y="Frequency")+
  style.print()
data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i)
sd(MArise)/(length(MArise)^0.5)
mean(MArise)
#3)
MArise <- (muttmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0))$deltadp
dataresrise <- contmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0) %>%
  group_by(x) %>%
  dplyr::summarize(CONrise=mean(deltadp))

source("~/Rfunction/style.print.R")
dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
  geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
  #scale_x_continuous(breaks = c(5,15,25))+
  labs(x="Average effect sizes of D-decreasing mutations",y="Frequency")+
  style.print()
data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i)
sd(MArise)/(length(MArise)^0.5)
mean(MArise)


library(plyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(Biostrings)
library(stringr)
library(tidyr)

# #1)calculate DP for all gene with one mutation (simulation)
load("~/project/HGT/exp/fig6.related/yeastgeneinf.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastcodonsupply.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
dpofref <- mclapply(mc.cores=6,1:nrow(infgene),function(x){
  gene <- infgene[x,1] 
  seq <- infgene[x,6]
  tmpyima <- which(substring(seq,seq(1,(nchar(seq)-2),by=3),seq(3,nchar(seq),by=3)) %in% c("TGA","TAG","TAA"))
  if(tmpyima[1] < nchar(seq)/3){
    myres1 <- data.frame(mutdpexp=99,mutdptai=99,datatype=99,wtDPexp=99,wtDPtai=99,gene,stringsAsFactors = F)
    
  } else {
    mutseq <- mclapply(mc.cores=10,4:(nchar(seq)-3),function(y){
      
      seqtmp <- seq
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
    mutseq$type[which(mutseq$seqall==seq)] <- "wt"
    alldp <- mclapply(mc.cores=10,1:nrow(mutseq),function(z){
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
          dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditai=(sum((yij-xijoftai)^2))^0.5) %>%
          dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptai=prod(Ditai)^(1/length(Ditai))) %>% 
          cbind(data.frame(stringsAsFactors = F,datatype))
        
      } else {
        res <- data.frame(stringsAsFactors = F,Dpexp=99, Dptai=99,datatype)
      }
      res
    }) %>% rbind.fill() %>% dplyr::filter(Dpexp!=99)
    myres1 <- (alldp %>% dplyr::filter(datatype=="mut"))
    myres0 <- alldp %>% dplyr::filter(datatype=="wt")
    myres1$wtDPexp <- myres0$Dpexp
    myres1$wtDPtai <- myres0$Dptai
    myres1$gene <- gene
    names(myres1)[1:2] <- c("mutdpexp","mutdptai")
  }
  myres1 
  
}) %>% rbind.fill() %>% dplyr::filter(wtDPexp != 99)
save(dpofref,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.ofrefgene.randommutation1.Rdata")


#2)calculate DP for MA genes
MAmutant1 <- read.csv("~/project/HGT/exp/model/MAmutant.csv",stringsAsFactors = F)
load("~/project/HGT/exp/fig6.related/yeastgeneinf.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastcodonsupply.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
filtergene <- mclapply(mc.cores=20,1:nrow(MAmutant1),function(x){
  atmp <- MAmutant1[x,]
  btmp <- infgene %>% dplyr::filter(genechr==atmp$chr & mintmp<atmp$pos & atmp$pos<maxtmp)
  if(nrow(btmp)>0){

    res <- cbind(btmp,atmp) %>% cbind(data.frame(type="T",x,stringsAsFactors = F))
  } else {

    res <- cbind(infgene[1,],atmp) %>% cbind(data.frame(type="F",x,stringsAsFactors = F))
  }
  res
}) %>% rbind.fill() %>% dplyr::filter(type=="T") %>%
  dplyr::filter(nchar(Reference)==nchar(Mutant))

dpofmut <- mclapply(mc.cores=20,1:nrow(filtergene),function(x){
  gene <- filtergene[x,1]
  exonreqion <- strsplit(filtergene[x,3],",")[[1]]
  exonsite1 <- lapply(1:length(exonreqion),function(y){
    a <- strsplit(exonreqion[y],"-")[[1]]
    if((min(as.numeric(a)) < filtergene[x,]$pos) & (max(as.numeric(a)) > filtergene[x,]$pos)){
      b <- data.frame(sta1=as.numeric(a)[1],end1=as.numeric(a)[2],len=abs(as.numeric(a)[1]-as.numeric(a)[2])+1,y,type="T",stringsAsFactors = F)
    } else {
      b <- data.frame(sta1=as.numeric(a)[1],end1=as.numeric(a)[2],len=abs(as.numeric(a)[1]-as.numeric(a)[2])+1,y,type="F",stringsAsFactors = F)
    }
    b
  }) %>% rbind.fill()

  if(exonsite1$sta1[1]-exonsite1$end1[1] < 0){
    exonsite <- exonsite1 %>% arrange(sta1)
    tmptype1 <- "zheng"
  } else{
    exonsite <- exonsite1 %>% arrange(desc(sta1))
    tmptype1 <- "fan"
  }

  if("T" %in% exonsite$type){
    at <- which(exonsite$type=="T")

    if(at==1){
      posoncds <- abs(exonsite$sta1[at]-filtergene[x,]$pos)+1
    } else {
      posoncds <- sum(exonsite$len[1:(at-1)]) + abs(exonsite$sta1[at]-filtergene[x,]$pos)+1
    }

    ref1 <- substr(filtergene[x,]$seq,posoncds,(posoncds+nchar(filtergene[x,]$Reference)-1))
    #filtergene[x,]
    if(ref1==filtergene[x,]$Reference){
      seq1 <- filtergene[x,]$seq
      str_sub(seq1,posoncds,(posoncds+nchar(filtergene[x,]$Reference)-1)) <- filtergene[x,]$Mutant
      tmptype2 <- "zheng"
    } else if(ref1==as.character(complement(DNAString(filtergene[x,]$Reference)))){
      seq1 <- filtergene[x,]$seq
      str_sub(seq1,posoncds,(posoncds+nchar(filtergene[x,]$Reference)-1)) <- as.character(complement(DNAString(filtergene[x,]$Mutant)))
      tmptype2 <- "fan"
    }
    #
    if(as.character(translate(DNAStringSet(seq1))[[1]]) == as.character(translate(DNAStringSet(filtergene[x,]$seq))[[1]])){
      muttype <- "Syn"
    } else {muttype <- "Non-Syn"}

    #calculate dp
    a <- substring(seq1,seq(1,(nchar(seq1)-2),by=3),seq(3,nchar(seq1),by=3)) %>% table() %>% as.data.frame()
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
        dplyr::summarize(Di.exp=(sum((yij-xij)^2)^0.5),Di.tai=(sum((yij-xijoftai)^2)^0.5)) %>%
        dplyr::summarize(DP.exp=round(prod(Di.exp)^(1/length(Di.exp)),3),DP.tai=round(prod(Di.tai)^(1/length(Di.tai)),3)) %>%
        cbind(data.frame(stringsAsFactors = F,gene))
      res1 <- res %>% cbind(data.frame(stringsAsFactors = F,tmptype1,tmptype2,muttype))
    } else {
      res1 <- data.frame(stringsAsFactors = F,DP.exp=99,DP.tai=99,gene,tmptype1=99,tmptype2=99,muttype=99)
    }
  } else {
    res1 <- data.frame(stringsAsFactors = F,DP.exp=99,DP.tai=99,gene,tmptype1=99,tmptype2=99,muttype=99)
  }
  res1
}) %>% rbind.fill() %>% dplyr::filter(DP.exp!=99)

save(dpofmut,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")

# # #2) plot 
##1)exp
load("/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.ofrefgene.randommutation1.Rdata")
dpofref <- dpofref[,c(1,4,6)]
names(dpofref)[1:2] <- c("mutdp","wtdp")
load("/mnt/data/home/chenfeng/project/HGT/exp/fig6.related/dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")
dpofmut <- dpofmut %>% dplyr::filter(muttype == "Syn")
dpofmut <- dpofmut[,c(1,3)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$gene,dpofref$gene)]

express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastmRNA.csv",stringsAsFactors = F) %>%
  arrange(desc(mRNA))
express <- express %>% dplyr::filter(ORF %in% dpofmut$gene) %>% group_by(ORF) %>%
  dplyr::filter(mRNA==max(mRNA)) %>% arrange(desc(mRNA))
##

i=0.1
experimentG <- express$ORF[1:floor(i*nrow(express))]
amutinf <- dpofmut %>% dplyr::filter(gene %in% experimentG) %>% group_by(gene) %>% dplyr::summarize(n=length(gene)) %>% group_by(n) %>% dplyr::summarize(nn=length(n))
amutinf1 <- dpofmut %>% dplyr::filter(gene %in% experimentG) %>% group_by(gene) %>% dplyr::summarize(n=length(gene))
mutif <- dpofmut %>% group_by(gene) %>% dplyr::summarize(nn=length(mutdp))

controlG <- lapply(1:length(amutinf$n), function(j){

  ainf <- amutinf[j,]
  ainf1 <- amutinf1 %>% dplyr::filter(n==ainf$n)
  binf <- mutif %>% dplyr::filter(nn==ainf$n & !(mutif$gene %in% experimentG))
  if(nrow(binf)==0){
    res <- data.frame(gene=ainf1$gene,nn=NA,Dp=NA)
  } else {
    binf$Dp <- dpofmut$wtdp[match(binf$gene,dpofmut$gene)]
    res <- (binf %>% arrange(Dp))[1:ainf$nn,]
  }
  res %>% cbind(type="Control")
}) %>% rbind.fill()

tmpgene <- data.frame(stringsAsFactors = F,type="Experiment",gene=setdiff(experimentG,controlG$gene)) %>%
  rbind(data.frame(stringsAsFactors = F,type="Control",gene=setdiff(controlG$gene,experimentG)))

k=1
typeoftwo <- unique(tmpgene$type)[k]
genetmp <- tmpgene$gene[which(tmpgene$type==typeoftwo)]

amut <- dpofmut %>% dplyr::filter(gene %in% genetmp)
bref <- dpofref %>% dplyr::filter(gene %in% genetmp)

atmp <- mclapply(mc.cores=30,1:length(genetmp), function(y){
  amutdata <- amut %>% dplyr::filter(gene==genetmp[y])
  brefall <- bref %>% dplyr::filter(gene==genetmp[y])

  subdata <- mclapply(mc.cores = 2,1:10000,function(x){
    set.seed(x)
    brefall[sample(nrow(brefall),nrow(amutdata),replace = T),] %>% cbind(data.frame(stringsAsFactors = F,x,type="con"))

  }) %>% rbind.fill()

  subdata %>% rbind(cbind(amutdata,data.frame(stringsAsFactors = F,x=10001,type="mut")))

}) %>% rbind.fill()

muttmp <- atmp %>% dplyr::filter(type=="mut")
contmp <- atmp %>% dplyr::filter(type=="con")

#1)
MArise <- (muttmp %>% dplyr::filter(mutdp>wtdp) %>% nrow())/nrow(muttmp)
dataresrise <- contmp %>%
  group_by(x) %>%
  dplyr::mutate(nn=length(x)) %>% as.data.frame() %>%
  group_by(x) %>%
  dplyr::summarize(CONrise=length(which(mutdp>wtdp))/unique(nn))

source("~/Rfunction/style.print.R")
dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
  geom_segment(aes(x=MArise,y=200,xend=MArise,yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
  scale_x_continuous(breaks = c(0.32,0.40,0.48,0.56,0.64))+
  labs(x="Percentage of D-increasing mutations \n among top 10% highly expressed genes",y="Frequency")+
  #labs(x="Percentage of D-increasing mutations \n among genes with 10% lowest D values",y="Frequency")+
  style.print()
data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > MArise)/10000,typeoftwo,i)
MArise
sd(dataresrise$CONrise)
mean(dataresrise$CONrise)
# #2)
MArise <- (muttmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0))$deltadp
dataresrise <- contmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0) %>%
  group_by(x) %>%
  dplyr::summarize(CONrise=mean(deltadp))

source("~/Rfunction/style.print.R")
dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
  geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
  scale_x_continuous(breaks = c(0.002,0.003,0.004,0.005))+
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
################
save(dpofmut,file = "~/project/HGT/exp/review1nd/data/dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")
save(dpofref,file = "~/project/HGT/exp/review1nd/data/dpexptai.ofrefgene.randommutation1.Rdata")

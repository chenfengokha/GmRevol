library(plyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(Biostrings)
library(stringr)
library(tidyr)

#1)calculate DP for all gene with one mutation (simulation)

load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
singlemutofrefseq <- mclapply(mc.cores=6,1:nrow(nameseq),function(x){
  myname <- nameseq$bname[x]
  myseq <- as.character(nameseq$myseq[x])
  mclapply(mc.cores=10,4:(nchar(myseq)-3),function(y){
    seqtmp <- myseq
    from <- substr(seqtmp,y,y)
    seqtmp1 <- as.character(translate(DNAString(seqtmp), no.init.codon=TRUE,if.fuzzy.codon="X"))
    to1 <- "A"
    str_sub(seqtmp,y,y) <- to1
    aseq <- seqtmp
    aseq1 <- as.character(translate(DNAString(aseq), no.init.codon=TRUE,if.fuzzy.codon="X"))

    to2 <- "T"
    str_sub(seqtmp,y,y) <- to2
    bseq <- seqtmp
    bseq1 <- as.character(translate(DNAString(bseq), no.init.codon=TRUE,if.fuzzy.codon="X"))

    to3 <- "G"
    str_sub(seqtmp,y,y) <- to3
    cseq <- seqtmp
    cseq1 <- as.character(translate(DNAString(cseq), no.init.codon=TRUE,if.fuzzy.codon="X"))

    to4 <- "C"
    str_sub(seqtmp,y,y) <- to4
    dseq <- seqtmp
    dseq1 <- as.character(translate(DNAString(dseq), no.init.codon=TRUE,if.fuzzy.codon="X"))
    restmp <- data.frame(stringsAsFactors = F,seqall=c(aseq,bseq,cseq,dseq),seqall1=c(aseq1,bseq1,cseq1,dseq1),from,to=c(to1,to2,to3,to4))
    restmp$type <- "wt"
    restmp$type[which(restmp$seqall1==seqtmp1 & restmp$seqall!=myseq)] <- "syn"
    restmp$type[which(restmp$seqall1!=seqtmp1)] <- "non"
    restmp[,-2]
  }) %>% rbind.fill() %>% unique() %>% dplyr::filter(type != "non") %>% cbind(data.frame(stringsAsFactors = F,myname))
}) %>% rbind.fill()

save(singlemutofrefseq,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.singlemutofrefseq.mg1655.1.Rdata")


load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.singlemutofrefseq.mg1655.1.Rdata")
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.Ecoli.Rdata")

mygene <- unique(singlemutofrefseq$myname)
dpofref <- mclapply(mc.cores=10,1:length(mygene),function(x){
  tmpa <- singlemutofrefseq %>% dplyr::filter(myname == mygene[x])
  
  alldp <- mclapply(mc.cores=3,1:nrow(tmpa),function(z){
    datatype <- tmpa[z,4]
    dataseq <- tmpa[z,1]
    myname <- tmpa[z,5]
    a <- substring(dataseq,seq(1,(nchar(dataseq)-2),by=3),seq(3,nchar(dataseq),by=3)) %>% table() %>% as.data.frame()
    names(a)[1] <- c("codon")
    a$codon <- as.vector(a$codon)
    a$amino <- codon$aa[match(a$codon,codon$codon)]
    
    xijdata <- xijofD
    xijdata$yfre <- a$Freq[match(xijdata$codon,a$codon)]
    xijdata$yfre[which(is.na(xijdata$yfre))] <- 0
    xijdata <- xijdata %>%
      group_by(amino) %>% dplyr::mutate(nn=sum(yfre)) %>% dplyr::mutate(yij=yfre/nn)
    
    xijdata$yij[which(is.nan(xijdata$yij))] <- 0
    xijdata %>% group_by(amino) %>%
      dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditaiexp304803=(sum((yij-xijoftai304803)^2))^0.5,Ditaiexp128812=(sum((yij-xijoftai128812)^2))^0.5,Ditaiexp100263=(sum((yij-xij100263)^2))^0.5) %>%
      dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptaiexp304803=prod(Ditaiexp304803)^(1/length(Ditaiexp304803)),Dptaiexp128812=prod(Ditaiexp128812)^(1/length(Ditaiexp128812)),Dptaiexp100263=prod(Ditaiexp100263)^(1/length(Ditaiexp100263))) %>%
      cbind(data.frame(stringsAsFactors = F,datatype,myname,from=tmpa[z,2],to=tmpa[z,3]))
    
  }) %>% rbind.fill()
  if(length(unique(alldp$datatype))>1){
    finres <- alldp %>% dplyr::filter(datatype=="syn") %>%
      cbind(data.frame(stringsAsFactors = F,WTdpexp=unique(alldp$Dpexp[which(alldp$datatype=="wt")]),
                       wtdptaiexp304803=unique(alldp$Dptaiexp304803[which(alldp$datatype=="wt")]),
                       wtdptaiexp128812=unique(alldp$Dptaiexp128812[which(alldp$datatype=="wt")]),
                       wtdptaiexp100263=unique(alldp$Dptaiexp100263[which(alldp$datatype=="wt")])))
    names(finres)[1:4] <- c("mutdpexp","mutdptaiexp304803","mutdptaiexp128812","mutdptaiexp100263")
  } else {
    finres1 <- data.frame(mutdpexp=99,datatype=99,WTdpexp=99,mutdptaiexp304803=99,mutdptaiexp128812=99,wtdptaiexp304803=99,wtdptaiexp128812=99,mutdptaiexp100263=99,wtdptaiexp100263=99,myname=99,from=99,to=99,stringsAsFactors = F)
  }
  
  finres
}) %>% rbind.fill() %>% dplyr::filter(mutdpexp!=99)

save(dpofref,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation.mg1655.Rdata")

#2)calculate DP for MA genes
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.Ecoli.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
MAmutant1 <- read.csv("~/project/HGT/exp/model/MG1655MA.csv",stringsAsFactors = F) %>% dplyr::filter(type %in% c("WT","MMR") & Amino.acid.mutation=="-" & bname %in% nameseq$bname)
dpofmut <- mclapply(mc.cores=10,1:nrow(MAmutant1),function(x){
  atmp <- MAmutant1[x,]
  myfrom <- substr(atmp$Codon.mutation,1,3)
  myto <- substr(atmp$Codon.mutation,6,8)
  #####
  basefromtmp <- strsplit(myfrom,"")[[1]]
  basetotmp <-  strsplit(myto,"")[[1]]
  nnrank <- which(basefromtmp != basetotmp)
  basefrom <- basefromtmp[nnrank]
  baseto <-  basetotmp[nnrank]
  ###
  btmp <- nameseq %>% dplyr::filter(toupper(bname)==toupper(atmp$bname)) %>% dplyr::mutate(len=nchar(myseq)) %>% dplyr::filter(len==max(len))
  dataseq <- btmp$myseq
  codonofgene <- substring(dataseq,seq(1,(nchar(dataseq)-2),by=3),seq(3,nchar(dataseq),by=3)) %>% table() %>% as.data.frame()
  names(codonofgene)[1] <- "codon"
  codonofgene$codon <- as.vector(codonofgene$codon)
  codonofgene %>% rbind(data.frame(stringsAsFactors = F,codon=c(myfrom,myto),Freq=c(-1,1))) -> codonofgene
  codonofgene$aa <- codon$aa[match(codonofgene$codon,codon$codon)]
  yfreMA <- codonofgene %>% dplyr::filter(!(aa %in% c(NA,"M","W"))) %>%
    group_by(aa,codon) %>% dplyr::summarize(newfre=sum(Freq))
  xijdata <- xijofD
  xijdata$yfre <- yfreMA$newfre[match(xijdata$codon,yfreMA$codon)]
  xijdata$yfre[which(is.na(xijdata$yfre))] <- 0
  xijdata %>%
    group_by(amino) %>% dplyr::mutate(aafre=sum(yfre)) %>%
    as.data.frame() %>% dplyr::mutate(yij=yfre/aafre) -> xijdata
  xijdata$yij[which(is.nan(xijdata$yij))] <- 0
  
  xijdata %>% group_by(amino) %>%
    dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditaiexp304803=(sum((yij-xijoftai304803)^2))^0.5,Ditaiexp128812=(sum((yij-xijoftai128812)^2))^0.5,Ditaiexp100263=(sum((yij-xij100263)^2))^0.5) %>%
    dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptaiexp304803=prod(Ditaiexp304803)^(1/length(Ditaiexp304803)),Dptaiexp128812=prod(Ditaiexp128812)^(1/length(Ditaiexp128812)),Dptaiexp100263=prod(Ditaiexp100263)^(1/length(Ditaiexp100263))) %>%
    cbind(data.frame(stringsAsFactors = F,gene=atmp$Gene,bname=atmp$bname,pos=atmp$Position,basefrom,baseto))
  
}) %>% rbind.fill()
save(dpofmut,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpofmutgene.MAline.MG1655.Rdata")
# #
# # #2)
##1)exp0.1 of D
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation.mg1655.Rdata")
dpofref <- dpofref[,c(1,6,9,7,8)]
names(dpofref)[c(1,3)] <- c("mutdp","wtdp")
dpofref$muttype <- paste(dpofref$from,dpofref$to)
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpofmutgene.MAline.MG1655.Rdata")
dpofmut <- dpofmut[,c(1,6,8,9)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$myname)]
dpofmut$muttype <- paste(dpofmut$basefrom,dpofmut$baseto)

##2)exptRNA of D
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation.mg1655.Rdata")
dpofref <- dpofref[,c(3,6,11,7,8)]
names(dpofref)[c(1,3)] <- c("mutdp","wtdp")
dpofref$muttype <- paste(dpofref$from,dpofref$to)
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpofmutgene.MAline.MG1655.Rdata")
dpofmut <- dpofmut[,c(3,6,8,9)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$myname)]
dpofmut$muttype <- paste(dpofmut$basefrom,dpofmut$baseto)

#########################################################################################
figtype <- c("endous","tRNA")
res <- lapply(1:2,function(xxx){
  figtmptype <- figtype[xxx]
  if(figtmptype == "endous"){
    load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation.mg1655.Rdata")
    dpofref <- dpofref[,c(1,6,9,7,8)]
    names(dpofref)[c(1,3)] <- c("mutdp","wtdp")
    dpofref$muttype <- paste(dpofref$from,dpofref$to)
    load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpofmutgene.MAline.MG1655.Rdata")
    dpofmut <- dpofmut[,c(1,6,8,9)]
    names(dpofmut)[1] <- "mutdp"
    dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$myname)]
    dpofmut$muttype <- paste(dpofmut$basefrom,dpofmut$baseto)
  } else {
    load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation.mg1655.Rdata")
    dpofref <- dpofref[,c(3,6,11,7,8)]
    names(dpofref)[c(1,3)] <- c("mutdp","wtdp")
    dpofref$muttype <- paste(dpofref$from,dpofref$to)
    load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpofmutgene.MAline.MG1655.Rdata")
    dpofmut <- dpofmut[,c(3,6,8,9)]
    names(dpofmut)[1] <- "mutdp"
    dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$myname)]
    dpofmut$muttype <- paste(dpofmut$basefrom,dpofmut$baseto)
  }
  
  express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
  express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2
  express <- express %>%
    arrange(desc(average))
  
  ranktmp = c(0.2,0.15,0.1,0.05,0.04,0.03)
  
  nuctmp <- data.frame(stringsAsFactors = F,s=c("A","G","C","T"),en=c("T","C","G","A"))
  
  resfin <- lapply(ranktmp,function(i){
    i=0.15
    experimentG <- express$log.TPM[1:floor(i*nrow(express))]
    
    amutinf <- dpofmut %>% dplyr::filter(bname %in% experimentG) %>% group_by(bname) %>% dplyr::summarize(n=length(bname)) %>% group_by(n) %>% dplyr::summarize(nn=length(n))
    
    mutif <- dpofmut %>% group_by(bname) %>% dplyr::summarize(nn=length(mutdp))
    
    controlG <- mclapply(mc.cores = 30,1:length(amutinf$n), function(j){
      
      ainf <- amutinf[j,]
      #ainf1 <- amutinf1 %>% dplyr::filter(n==ainf$n)
      binf <- mutif %>% dplyr::filter(nn==ainf$n & !(mutif$bname %in% experimentG))
      if(nrow(binf)==0){
        res <- data.frame(bname=NA,nn=NA,Dp=NA)
      } else {
        binf$Dp <- dpofmut$wtdp[match(binf$bname,dpofmut$bname)]
        res <- (binf %>% arrange((Dp)))[1:ainf$nn,]
      }
      res %>% cbind(type="Control")
    }) %>% rbind.fill() %>% dplyr::filter(!is.na(Dp))
    
    tmpgene <- data.frame(stringsAsFactors = F,type="Experiment",gene=setdiff(experimentG,controlG$bname)) %>%
      rbind(data.frame(stringsAsFactors = F,type="Control",gene=setdiff(controlG$bname,experimentG)))
    
    mclapply(mc.cores = 1,1:2, function(k){
      
      typeoftwo <- unique(tmpgene$type)[k]
      genetmp <- tmpgene$gene[which(tmpgene$type==typeoftwo)]
      
      
      amut <- dpofmut %>% dplyr::filter(bname %in% genetmp)
      bref <- dpofref %>% dplyr::filter(myname %in% genetmp)
      
      atmp <- mclapply(mc.cores=30,1:length(unique(amut$bname)), function(y){
        amutdata <- amut %>% dplyr::filter(bname==unique(amut$bname)[y])
        brefall <- bref %>% dplyr::filter(myname==unique(amut$bname)[y])
        names(brefall)[2] <- "bname"
        
        subdata <- lapply(1:nrow(amutdata),function(x){
          tmpdata <- brefall %>% dplyr::filter((muttype %in% c(amutdata$muttype[x],paste(nuctmp$en[which(nuctmp$s == amutdata$myfrom[x])],nuctmp$en[which(nuctmp$s == amutdata$myto[x])]))))
          
          mclapply(mc.cores = 3,1:10000,function(z){
            set.seed(z)
            tmpdata[sample(nrow(tmpdata),1),][c(-5,-4)] %>% cbind(data.frame(stringsAsFactors = F,z,type="con"))
            
          }) %>% rbind.fill()
        }) %>% rbind.fill()
        
        subdata %>% rbind(cbind(amutdata[c(-3,-4)],data.frame(stringsAsFactors = F,z=10001,type="mut")))
        
      }) %>% rbind.fill()
      
      muttmp <- atmp %>% dplyr::filter(type=="mut")
      contmp <- atmp %>% dplyr::filter(type=="con")
      ######################res
      #1)
      MArise <- (muttmp %>% dplyr::filter(mutdp>wtdp)  %>% nrow())/nrow(muttmp)
      dataresrise <- contmp %>%
        group_by(z) %>%
        dplyr::mutate(nn=length(z)) %>% as.data.frame() %>%
        group_by(z) %>%
        dplyr::summarize(CONrise=length(which(mutdp>wtdp))/unique(nn))
      
      source("~/Rfunction/style.print.R")
      dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
        geom_segment(aes(x=MArise,y=200,xend=MArise,yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
        #scale_x_continuous(breaks = c(0.35,0.45,0.55))+
        #labs(x="Percentage of D-increasing mutations \namong top 15% highly expressed genes",y="Frequency")+
        labs(x="Percentage of D−increasing mutations\namong genes with 15% lowest D values",y="Frequency")+
        style.print()
      write.csv(dataresrise, file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig. S10bL.csv")
      data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > MArise)/10000,typeoftwo,i,fig="a")
      # typeoftwo
      
      MArise
      sd(dataresrise$CONrise)
      mean(dataresrise$CONrise)
      
      #2)
      MArise <- (muttmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0))$deltadp
      dataresrise <- contmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0) %>%
        group_by(z) %>%
        dplyr::summarize(CONrise=mean(deltadp))
      
      source("~/Rfunction/style.print.R")
      dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
        geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
        #scale_x_continuous(limits = c(0.002,0.005))+
        labs(x="Average effect sizes of D-increasing mutations",y="Frequency")+
        style.print()
      write.csv(dataresrise, file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig. S12b.csv")
      data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i,fig="b")
      sd(MArise)/(length(MArise)^0.5)
      mean(MArise)
      #3)
      MArise <- (muttmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0))$deltadp
      dataresrise <- contmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0) %>%
        group_by(z) %>%
        dplyr::summarize(CONrise=mean(deltadp))
      
      source("~/Rfunction/style.print.R")
      dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
        geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
        #scale_x_continuous(breaks = c(5,15,25))+
        labs(x="Average effect sizes of D-decreasing mutations",y="Frequency")+
        style.print()
      write.csv(dataresrise, file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig. S12c.csv")
      data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i,fig="c")
      sd(MArise)/(length(MArise)^0.5)
      mean(MArise)
      rbind(res1,res2,res3) %>% cbind(data.frame(stringsAsFactors = F,n=length(unique(amut$bname))))
      
    }) %>% rbind.fill()
    
  }) %>% rbind.fill()
  resfin$figtmp <- figtmptype
  resfin
}) %>% rbind.fill()

res %>% dplyr::filter(typeoftwo == "Experiment") %>% dplyr::filter(fig == "a")









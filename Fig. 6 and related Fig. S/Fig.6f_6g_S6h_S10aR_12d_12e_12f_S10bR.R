library(plyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(Biostrings)
library(stringr)
library(tidyr)

# #1)calculate DP for all gene with one mutation (simulation)
load("~/project/HGT/exp/fig6.related/yeastgeneinf.Rdata")

singlemutofrefseq <- mclapply(mc.cores=6,1:nrow(infgene),function(x){
  gene <- infgene[x,1]
  seq <- infgene[x,6]
  tmpyima <- which(substring(seq,seq(1,(nchar(seq)-2),by=3),seq(3,nchar(seq),by=3)) %in% c("TGA","TAG","TAA"))
  if(tmpyima[1] < nchar(seq)/3){
    mutseq <- data.frame(seqall=99,type=99,gene,stringsAsFactors = F)

  } else {

    mutseq <- mclapply(mc.cores=10,4:(nchar(seq)-3),function(y){
      seqtmp <- seq
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
      restmp$type[which(restmp$seqall1==seqtmp1 & restmp$seqall!=seq)] <- "syn"
      restmp$type[which(restmp$seqall1!=seqtmp1)] <- "non"
      restmp[,-2]
    }) %>% rbind.fill() %>% unique() %>% dplyr::filter(type != "non") %>% cbind(data.frame(stringsAsFactors = F,gene))
  }
  mutseq
}) %>% rbind.fill()
save(singlemutofrefseq,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.singlemutofrefseq.yeast.Rdata")

# #length(mygene)
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.singlemutofrefseq.yeast.Rdata")
singlemutofrefseq <- singlemutofrefseq %>% dplyr::filter(type!="99")
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.SC.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
mygene <- unique(singlemutofrefseq$gene)

dpofref <- mclapply(mc.cores=15,1:length(mygene),function(x){
  tmpa <- singlemutofrefseq %>% dplyr::filter(gene == mygene[x])
  
  alldp <- mclapply(mc.cores=6,1:nrow(tmpa),function(z){
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
      dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditaiexp243449=(sum((yij-xij243449)^2))^0.5,Ditaiexp21558172=(sum((yij-xij21558172)^2))^0.5,Ditaiexp33581077=(sum((yij-xij33581077)^2))^0.5) %>%
      dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptaiexp243449=prod(Ditaiexp243449)^(1/length(Ditaiexp243449)),Dptaiexp21558172=prod(Ditaiexp21558172)^(1/length(Ditaiexp21558172)),Dptaiexp33581077=prod(Ditaiexp33581077)^(1/length(Ditaiexp33581077))) %>%
      cbind(data.frame(stringsAsFactors = F,datatype,from=tmpa[z,2],to=tmpa[z,3]))
    
  }) %>% rbind.fill()
  if(length(unique(alldp$datatype))>1){
    finres <- alldp %>% dplyr::filter(datatype=="syn") %>%
      cbind(data.frame(stringsAsFactors = F,WTdpexp=unique(alldp$Dpexp[which(alldp$datatype=="wt")]),wtdptaiexp243449=unique(alldp$Dptaiexp243449[which(alldp$datatype=="wt")]),
                       wtdptaiexp21558172=unique(alldp$Dptaiexp21558172[which(alldp$datatype=="wt")]),wtdptaiexp33581077=unique(alldp$Dptaiexp33581077[which(alldp$datatype=="wt")])))
    
    names(finres)[1:4] <- c("mutdpexp","mutdptaiexp243449","mutdptaiexp21558172","mutdptaiexp33581077")
  } else {
    finres <- data.frame(mutdpexp=99,mutdptaiexp243449=99,mutdptaiexp21558172=99,mutdptaiexp33581077=99,datatype=99,WTdpexp=99,wtdptaiexp243449=99,wtdptaiexp21558172=99,wtdptaiexp33581077=99,from=99,to=99,stringsAsFactors = F)
  }
  
  finres$gene <- mygene[x]
  
  finres
  
}) %>% rbind.fill() %>% dplyr::filter(WTdpexp != 99)
save(dpofref,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation1.Rdata")



# #2)calculate DP for MA genes
MAmutant1 <- read.csv("~/project/HGT/exp/model/MAmutant.csv",stringsAsFactors = F)
load("~/project/HGT/exp/fig6.related/yeastgeneinf.Rdata")
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.SC.Rdata")
codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
filtergene <- mclapply(mc.cores=50,1:nrow(MAmutant1),function(x){
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
  myfrom <- filtergene$Reference[x]
  myto <- filtergene$Mutant[x]
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
    
    if(ref1==filtergene[x,]$Reference){
      seq1 <- filtergene[x,]$seq
      str_sub(seq1,posoncds,(posoncds+nchar(filtergene[x,]$Reference)-1)) <- filtergene[x,]$Mutant
      tmptype2 <- "zheng"
    } else if(ref1==as.character(Biostrings::complement(DNAString(filtergene[x,]$Reference)))){
      seq1 <- filtergene[x,]$seq
      str_sub(seq1,posoncds,(posoncds+nchar(filtergene[x,]$Reference)-1)) <- as.character(Biostrings::complement(DNAString(filtergene[x,]$Mutant)))
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
      xijdata <- xijofD
      xijdata$yfre <- a$Freq[match(xijdata$codon,a$codon)]
      xijdata$yfre[which(is.na(xijdata$yfre))] <- 0
      xijdata <- xijdata %>%
        group_by(amino) %>% dplyr::mutate(nn=sum(yfre)) %>% dplyr::mutate(yij=yfre/nn)
      
      
      xijdata$yij[which(is.nan(xijdata$yij))] <- 0
      res <- xijdata %>% group_by(amino) %>%
        dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditaiexp243449=(sum((yij-xij243449)^2))^0.5,Ditaiexp21558172=(sum((yij-xij21558172)^2))^0.5,Ditaiexp33581077=(sum((yij-xij33581077)^2))^0.5) %>%
        dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptaiexp243449=prod(Ditaiexp243449)^(1/length(Ditaiexp243449)),Dptaiexp21558172=prod(Ditaiexp21558172)^(1/length(Ditaiexp21558172)),Dptaiexp33581077=prod(Ditaiexp33581077)^(1/length(Ditaiexp33581077))) %>%
        cbind(data.frame(stringsAsFactors = F,gene))
      res1 <- res %>% cbind(data.frame(stringsAsFactors = F,tmptype1,tmptype2,muttype))
    } else {
      res1 <- data.frame(stringsAsFactors = F,Dpexp=99,Dptaiexp243449=99,Dptaiexp21558172=99,Dptaiexp33581077=99,gene=99,tmptype1=99,tmptype2=99,muttype=99)
    }
  } else {
    res1 <- data.frame(stringsAsFactors = F,Dpexp=99,Dptaiexp243449=99,Dptaiexp21558172=99,Dptaiexp33581077=99,gene=99,tmptype1=99,tmptype2=99,muttype=99)
  }
  res1 %>% cbind(data.frame(stringsAsFactors = F,myfrom,myto))
}) %>% rbind.fill() %>% dplyr::filter(Dpexp!=99)

save(dpofmut,file = "/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")

# # # #2) plot
##1)exp0.1 of D
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation1.Rdata")
dpofref <- dpofref[,c(1,8,12,6,7)]
names(dpofref)[c(1,2)] <- c("mutdp","wtdp")
dpofref$muttype <- paste(dpofref$from,dpofref$to)
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")
dpofmut <- (dpofmut %>% dplyr::filter(muttype == "Syn"))[,c(1,5,9,10)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$gene,dpofref$gene)]
dpofmut$muttype <- paste(dpofmut$myfrom,dpofmut$myto)


##2)exptRNA of D
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.ofrefgene.randommutation1.Rdata")
dpofref <- dpofref[,c(3,10,12,6,7)]
names(dpofref)[c(1,2)] <- c("mutdp","wtdp")
dpofref$muttype <- paste(dpofref$from,dpofref$to)
load("/mnt/data/home/chenfeng/project/chenfengdata5/NC2025review2nd/tRNAsupply/fig6.dpexptai.mutgene.MAline.3paper.allcodon.1.Rdata")
dpofmut <- (dpofmut %>% dplyr::filter(muttype == "Syn"))[,c(3,5,9,10)]
names(dpofmut)[1] <- "mutdp"
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$gene,dpofref$gene)]
dpofmut$muttype <- paste(dpofmut$myfrom,dpofmut$myto)


express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastmRNA.csv",stringsAsFactors = F) %>%
  arrange(desc(mRNA))

ranktmp = c(0.2,0.15,0.1,0.05,0.04,0.03)

nuctmp <- data.frame(stringsAsFactors = F,s=c("A","G","C","T"),en=c("T","C","G","A"))
####################################
resfin <- lapply(ranktmp,function(i){
  i=0.15
  experimentG <- express$ORF[1:floor(i*nrow(express))]
  
  amutinf <- dpofmut %>% dplyr::filter(gene %in% experimentG) %>% group_by(gene) %>% dplyr::summarize(n=length(gene)) %>% group_by(n) %>% dplyr::summarize(nn=length(n))
  amutinf1 <- dpofmut %>% dplyr::filter(gene %in% experimentG) %>% group_by(gene) %>% dplyr::summarize(n=length(gene))
  mutif <- dpofmut %>% group_by(gene) %>% dplyr::summarize(nn=length(mutdp))
  
  controlG <- mclapply(mc.cores = 30,1:length(amutinf$n), function(j){
    
    ainf <- amutinf[j,]
    ainf1 <- amutinf1 %>% dplyr::filter(n==ainf$n)
    binf <- mutif %>% dplyr::filter(nn==ainf$n & !(mutif$gene %in% experimentG))
    if(nrow(binf)==0){
      res <- data.frame(gene=ainf1$gene,nn=NA,Dp=NA)
    } else {
      binf$Dp <- dpofmut$wtdp[match(binf$gene,dpofmut$gene)]
      res <- (binf %>% arrange((Dp)))[1:ainf$nn,]
    }
    res %>% cbind(type="Control")
  }) %>% rbind.fill()
  
  tmpgene <- data.frame(stringsAsFactors = F,type="Experiment",gene=setdiff(experimentG,controlG$gene)) %>%
    rbind(data.frame(stringsAsFactors = F,type="Control",gene=setdiff(controlG$gene,experimentG)))
  
  mclapply(mc.cores = 2,1:2, function(k){
    k=2
    typeoftwo <- unique(tmpgene$type)[k]
    genetmp <- tmpgene$gene[which(tmpgene$type==typeoftwo)]
    
    amut <- dpofmut %>% dplyr::filter(gene %in% genetmp)
    bref <- dpofref %>% dplyr::filter(gene %in% genetmp)
    
    atmp <- mclapply(mc.cores=30,1:length(unique(amut$gene)), function(y){
      amutdata <- amut %>% dplyr::filter(gene==unique(amut$gene)[y])
      brefall <- bref %>% dplyr::filter(gene==unique(amut$gene)[y])
      
      subdata <- lapply(1:nrow(amutdata),function(x){
        
        tmpdata <- brefall %>% dplyr::filter((muttype %in% c(amutdata$muttype[x],paste(nuctmp$en[which(nuctmp$s == amutdata$myfrom[x])],nuctmp$en[which(nuctmp$s == amutdata$myto[x])]))))
        if(nrow(tmpdata)  > 1){
          tmpres <- mclapply(mc.cores = 2,1:10000,function(z){
            set.seed(z)
            tmpdata[sample(nrow(tmpdata),1),][c(-5,-4)] %>% cbind(data.frame(stringsAsFactors = F,z,type="con"))
            
          }) %>% rbind.fill()
        } else {tmpres <- data.frame(stringsAsFactors = F,mutdp=99,wtdp=99,gene=99,muttype=99,z=99,type)}
        tmpres
      }) %>% rbind.fill()
      
      subdata %>% rbind(cbind(amutdata[,c(-3,-4)],data.frame(stringsAsFactors = F,z=10001,type="mut")))
      
    }) %>% rbind.fill() %>% dplyr::filter(mutdp!=99)
    
    
    muttmp <- atmp %>% dplyr::filter(type=="mut")
    contmp <- atmp %>% dplyr::filter(type=="con")
    
    #1)
    MArise <- (muttmp %>% dplyr::filter(mutdp>wtdp) %>% nrow())/nrow(muttmp)
    dataresrise <- contmp %>%
      group_by(z) %>%
      dplyr::mutate(nn=length(z)) %>% as.data.frame() %>%
      group_by(z) %>%
      dplyr::summarize(CONrise=length(which(mutdp>wtdp))/unique(nn))
    
    source("~/Rfunction/style.print.R")
    dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
      geom_segment(aes(x=MArise,y=200,xend=MArise,yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
      scale_x_continuous(breaks = c(0.32,0.40,0.48,0.56,0.64))+
      labs(x="Percentage of D-increasing mutations \n among top 15% highly expressed genes",y="Frequency")+
      #labs(x="Percentage of D-increasing mutations \n among genes with 15% lowest D values",y="Frequency")+
      style.print()
    write.csv(dataresrise,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S10bR.csv")
    data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > MArise)/10000,typeoftwo,i,fig="a")
    MArise
    sd(dataresrise$CONrise)
    mean(dataresrise$CONrise)
    # #2)
    MArise <- (muttmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0))$deltadp
    dataresrise <- contmp %>% dplyr::mutate(deltadp=mutdp-wtdp) %>% dplyr::filter(deltadp>0) %>%
      group_by(z) %>%
      dplyr::summarize(CONrise=mean(deltadp))
    
    source("~/Rfunction/style.print.R")
    dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
      geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
      #scale_x_continuous(breaks = c(0.002,0.003,0.004,0.005))+
      labs(x="Average effect sizes of D-increasing mutations",y="Frequency")+
      style.print()
    write.csv(dataresrise,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S12e.csv")
    data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i,fig="b")
    sd(MArise)/(length(MArise)^0.5)
    mean(MArise)
    #3)
    MArise <- (muttmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0))$deltadp
    dataresrise <- contmp %>% dplyr::mutate(deltadp=wtdp-mutdp) %>% dplyr::filter(deltadp>0) %>%
      group_by(z) %>%
      dplyr::summarize(CONrise=mean(deltadp))
    # 
    source("~/Rfunction/style.print.R")
    dataresrise %>% ggplot(aes(x=CONrise))+geom_histogram(bins = 30)+
      geom_segment(aes(x=mean(MArise),y=200,xend=mean(MArise),yend=20),arrow=arrow(length = unit(0.5,"cm")),color="red")+
      #scale_x_continuous(breaks = c(5,15,25))+
      labs(x="Average effect sizes of D-decreasing mutations",y="Frequency")+
      style.print()
    write.csv(dataresrise,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S12f.csv")
    data.frame(stringsAsFactors = F,P=sum(dataresrise$CONrise > mean(MArise))/10000,typeoftwo,i,fig="c")
    sd(MArise)/(length(MArise)^0.5)
    mean(MArise)
    
    rbind(res1,res2,res3) %>% cbind(data.frame(stringsAsFactors = F,n=length(experimentG)))
    
  }) %>% rbind.fill()
  
}) %>% rbind.fill()
resfin



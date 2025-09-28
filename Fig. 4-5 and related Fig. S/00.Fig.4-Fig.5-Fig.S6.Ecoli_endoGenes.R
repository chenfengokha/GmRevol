library(ppcor);
library(scales);
library(Biostrings)
##calculat DP
ecolicds <- readDNAStringSet("~/project/HGT/exp/paper/EcoliCDS.fasta")
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2
top.1exp <- sort(express$average,decreasing = T)[floor(0.1*nrow(express))]

nameseq <- mclapply(mc.cores=20,1:length(ecolicds),function(x){
  myname1 <- strsplit(names(ecolicds)[x]," ")[[1]][2]
  myname <- substr(myname1,7,(nchar(myname1)-1))
  bname1 <- strsplit(names(ecolicds)[x]," ")[[1]][3]
  bname <- (substr(bname1,12,16))
  myseq <- as.character(ecolicds[[x]])
  data.frame(stringsAsFactors = F,myname,myseq,bname,len=nchar(myseq))

}) %>% rbind.fill() %>% group_by(bname) %>% dplyr::filter(len==max(len) & bname %in% express$log.TPM) %>% unique()

codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
sequencedata <- mclapply(mc.cores=20,1:nrow(nameseq),function(x){
  seq <- substring(as.character(nameseq$myseq[x]),seq(1,(nchar(nameseq$myseq[x])-2),by=3),seq(3,nchar(nameseq$myseq[x]),by=3))
  a <- as.data.frame(table(seq));
  names(a) <- c("codon","fre");
  a$codon <- as.vector(a$codon)
  a$gene <- nameseq$bname[x]
  a$amino <- codon$aa[match(a$codon,codon$codon)]
  a;
}) %>% rbind.fill()

sequencedata$exp <- express$average[match(sequencedata$gene,express$log.TPM)]

xijofD <- sequencedata %>% 
  dplyr::filter(exp>top.1exp) %>%
  group_by(codon,amino) %>% dplyr::summarise(freq=sum(fre*exp)) %>% as.data.frame() %>%
  dplyr::filter(amino != "*" & amino != "W" & amino != "M") %>%
  group_by(amino) %>% dplyr::mutate(nn=sum(freq)) %>% dplyr::mutate(xij=freq/nn)
load("/home/chenfeng/project/HGT/exp/tAIcodecoli.Rdata")
xijofD$tai <- tAIofecoli$wi[match(xijofD$codon,tAIofecoli$codon)]

xijofD %>% group_by(amino) %>% dplyr::mutate(nntai=sum(tai)) %>% dplyr::mutate(xijoftai=tai/nntai) -> xijofD
xijofD %>% group_by(amino) %>% dplyr::summarize(nxij=sum(xij),nxijoftai=sum(xijoftai)) %>% as.data.frame()

DP.exp.tai <- mclapply(mc.cores = 10,1:nrow(nameseq),function(x){
  nucleotidedf <- strsplit(as.character(nameseq$myseq[x]),"")[[1]]
  mydf0 <- substring(as.character(nameseq$myseq[x]),seq(1,(nchar(nameseq$myseq[x])-2),by=3),seq(3,nchar(nameseq$myseq[x]),by=3))
  mydf <- mydf0 %>% table() %>% as.data.frame()
  names(mydf)[1] <- c("codon")
  mydf$codon <- as.vector(mydf$codon)
  xij <- xijofD
  xij$yfre <- mydf$Freq[match(xij$codon,mydf$codon)]
  xij$yfre[which(is.na(xij$yfre))] <- 0
  xij %>% group_by(amino) %>% dplyr::mutate(ny=sum(yfre)) %>% group_by(amino,codon) %>% dplyr::mutate(yij=yfre/ny) -> xij
  xij$yij[which(xij$yij=="NaN")] <- 0
  xij %>%
    group_by(amino) %>% 
    dplyr::summarize(Di.exp=(sum((yij-xij)^2)^0.5),Di.tai=(sum((yij-xijoftai)^2)^0.5)) %>% 
    dplyr::summarize(DP.exp=round(prod(Di.exp)^(1/length(Di.exp)),3),DP.tai=round(prod(Di.tai)^(1/length(Di.tai)),3)) %>% cbind(nameseq[x,c(1,3:4)]) %>%
    cbind(data.frame(stringsAsFactors = F,Apercentage=length(which(nucleotidedf=="A"))/length(nucleotidedf),
                     Tpercentage=length(which(nucleotidedf=="T"))/length(nucleotidedf),
                     Gpercentage=length(which(nucleotidedf=="G"))/length(nucleotidedf),
                     Cpercentage=length(which(nucleotidedf=="C"))/length(nucleotidedf),
                     GCpercentage=length(which(nucleotidedf=="G" | nucleotidedf=="C"))/length(nucleotidedf)))
  
}) %>% rbind.fill() %>% unique()

##gene importance
genefitness <- read.csv("~/project/HGT/exp/paper/genefitness.PMID32994326.csv",header = TRUE,stringsAsFactors = F)
alldata <- mclapply(mc.cores=20,3:ncol(genefitness),function(x){
  a <- genefitness[,c(1,2,x)]
  names(a)[3] <- "growth"
  a$importance <- 1-a$growth
  a
}) %>% rbind.fill() %>%
  group_by(B.numBers,Gene) %>%
  dplyr::summarize(imp=median(importance),n=length(importance)) %>%
  as.data.frame()

DP.exp.tai$importance <- alldata$imp[match(toupper(DP.exp.tai$bname),alldata$B.numBers)]
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2


DP.exp.tai$mRNAexp <- express$average[match(DP.exp.tai$bname,express$log.TPM)]
save(DP.exp.tai,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpexptai.ecoligene.Rdata")

#################################################################################################################################################################
#Fig 4
library(tidyr)
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpexptai.ecoligene.Rdata")
rnaThresQuant <- c(1,0.5,0.3,0.1,0.05,0.04,0.03)

geneDP0 <- DP.exp.tai %>% arrange(mRNAexp) %>% dplyr::filter(!if_any(c("mRNAexp","importance"), is.na))
names(geneDP0)[1] <- "DP"
dfcor.D.imp <- lapply(rnaThresQuant,function(x){
  thisSubset1 <- (geneDP0 %>% arrange(desc(mRNAexp)))[1:floor(x*nrow(geneDP0)),] %>% dplyr::filter(!if_any(c("importance","mRNAexp"), is.na))
  # # plot
  # thisSubset1 %>%
  #   ggplot(aes(x=DP,y=importance))+
  #   geom_point(shape=1,size=0.2)+
  #   scale_y_continuous(trans = "log2",breaks = c(0.002,0.02,0.2))+
  #   geom_smooth(method = "lm",se=F)+
  #   labs(x="D",y="Importance")+style.print()
  # cor.test(thisSubset1$DP,thisSubset1$importance,alternative="less",method = "s")
  res <- lapply(1:1000, function(y){
    set.seed(y)
    thisSubset <- thisSubset1[sample(1:nrow(thisSubset1),nrow(thisSubset1),replace = T),]
    set.seed(y)
    randomSubset <- geneDP0[sample(1:nrow(geneDP0),nrow(thisSubset1),replace = T),]  %>% dplyr::filter(!if_any(c("importance","mRNAexp"), is.na))
    pcorObj1 <- cor.test(thisSubset$DP,thisSubset$importance,method="s",alternative="less");
    pcorObj2 <- cor.test(randomSubset$DP,randomSubset$importance,method="s",alternative="less");
    data.frame(stringsAsFactors = F,
               thres = rep(paste(x*100,"%",sep = ""),2),
               estimate = c(pcorObj1$estimate,pcorObj2$estimate),
               p.val = c(pcorObj1$p.value,pcorObj2$p.value),
               
               corWith = c("real","random"))
  }) %>% rbind.fill()
  res %>% group_by(thres,corWith) %>% dplyr::summarize(mrho=mean(estimate),sd=sd(estimate)) %>% 
    cbind(data.frame(stringsAsFactors = F,pvalue=rep(wilcox.test((res %>% dplyr::filter(corWith == "real"))$estimate,(res %>% dplyr::filter(corWith == "random"))$estimate)$p.value,2)))
  
}) %>% rbind.fill()

source("~/Rfunction/style.print.R")
dfcor.D.imp$thres <- factor(dfcor.D.imp$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.D.imp$corWith <- factor(dfcor.D.imp$corWith,levels = c("real","random"))
dfcor.D.imp %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (Importance ~ D)",limits = c(-0.4,0.1),breaks = c(-0.4,-0.2,0,0.2,0.4,0.6,0.8)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes (top expression percentile)") +
  style.print()

##########################################################################################
#### fig 5
library(ggplot2)
library(ggExtra)

load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpexptai.ecoligene.Rdata")

rnaThresQuant <- c(1,0.5,0.3,0.1,0.05,0.04,0.03)
geneDP0 <- DP.exp.tai %>% arrange(mRNAexp) %>% dplyr::filter(!if_any(c("mRNAexp"), is.na))
names(geneDP0)[1] <- "DP"

##partial cor of A%,T%,G%,C%,and GC%
# names(geneDP0)[c(9,10)] <- c("d","Apercentage")
# head(geneDP0)
dfcor.exp.D <- mclapply(mc.cores = 7,rnaThresQuant,function(x){
  thisSubset1 <- (geneDP0 %>% arrange(desc(mRNAexp)))[1:floor(x*nrow(geneDP0)),] %>% dplyr::filter(!if_any(c("mRNAexp"), is.na))
  # # #plot
  # p1 <- thisSubset1 %>%
  #   ggplot(aes(x=DP,y=mRNAexp))+
  #   geom_point(shape=1,size=0.2)+
  #   scale_y_continuous(breaks = c(9,11,13))+
  #   geom_smooth(method = "lm",se=F)+
  #   labs(x="D",y="Expression")+style.print()
  # ggMarginal(p1, margins = "y", type = "histogram", groupColour = F, groupFill = F)
  # cor.test(thisSubset1$DP,thisSubset1$mRNAexp,alternative="greater",method = "s")
  #
  res <- lapply(1:1000, function(y){
    set.seed(y)
    thisSubset <- thisSubset1[sample(1:nrow(thisSubset1),nrow(thisSubset1),replace = T),]
    set.seed(y)
    randomSubset <- geneDP0[sample(1:nrow(geneDP0),nrow(thisSubset1),replace = T),] %>% dplyr::filter(!if_any(c("mRNAexp"), is.na))
    pcorObj1 <- cor.test(thisSubset$DP,thisSubset$mRNAexp,method="s",alternative="greater");
    pcorObj2 <- cor.test(randomSubset$DP,randomSubset$mRNAexp,method="s",alternative="greater");
    pcorA1 <- pcor.test(thisSubset$mRNAexp,thisSubset$DP,thisSubset$Apercentage,method="s")
    pcorA2 <- pcor.test(randomSubset$mRNAexp,randomSubset$DP,randomSubset$Apercentage,method="s")
    data.frame(stringsAsFactors = F,
               thres = rep(paste(x*100,"%",sep = ""),2),
               estimate = c(pcorObj1$estimate,pcorObj2$estimate),
               p.val = c(pcorObj1$p.value,pcorObj2$p.value),
               Aestimate = c(pcorA1$estimate,pcorA2$estimate),
               Ap.val = c(pcorA1$p.value,pcorA2$p.value),
               corWith = c("real","random"))
  }) %>% rbind.fill()
  res %>% group_by(thres,corWith) %>% dplyr::summarize(mrho=mean(estimate),sd=sd(estimate)) %>% 
    cbind(data.frame(stringsAsFactors = F,pvalue=rep(wilcox.test((res %>% dplyr::filter(corWith == "real"))$estimate,(res %>% dplyr::filter(corWith == "random"))$estimate)$p.value,2)))
}) %>% rbind.fill()

source("~/Rfunction/style.print.R")
dfcor.exp.D$thres <- factor(dfcor.exp.D$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.exp.D$corWith <- factor(dfcor.exp.D$corWith,levels = c("real","random"))
dfcor.exp.D %>%
  ggplot(aes(y=mrho,x=thres,fill=corWith,group=corWith)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8) +
  geom_errorbar(aes(ymax=mrho+sd,ymin=mrho-sd),position = position_dodge(width = 0.8),width = 0.3)+
  scale_fill_manual(values=c("#9933FF", "grey"))+
  scale_y_continuous("rho (Expression ~ D)",limits = c(-0.25,0.6),breaks = c(-0.4,-0.2,0,0.2,0.4,0.6,0.8)) +
  #scale_y_continuous("Partial rho(Expression ~ D) controlling for\nbase composition of nucleotide GC",limits = c(-0.4,0.6),breaks = c(-0.4,-0.2,0,0.2,0.4,0.6,0.8)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes\n(top expression percentile)") +
  style.print()

library(ppcor);
library(scales);
###Ecoli as an example
######load refinf
ecolicds <- readDNAStringSet("~/project/HGT/exp/paper/EcoliCDS.fasta")
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2
express <- express %>% arrange(desc(average))
#genename and cds transform
nameseq <- mclapply(mc.cores=20,1:length(ecolicds),function(x){
  myname1 <- strsplit(names(ecolicds)[x]," ")[[1]][2]
  myname <- substr(myname1,7,(nchar(myname1)-1))
  bname1 <- strsplit(names(ecolicds)[x]," ")[[1]][3]
  bname <- (substr(bname1,12,16))
  myseq <- as.character(ecolicds[[x]])
  data.frame(stringsAsFactors = F,myname,myseq,bname,len=nchar(myseq))
  
}) %>% rbind.fill() %>% group_by(bname) %>% dplyr::filter(len==max(len) & bname %in% express$log.TPM) %>% unique()
save(nameseq,file = "/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")


codon <- read.table("/mnt/data/home/chenfeng/project/codonpaper/expriment/fcs/codon.txt",header = TRUE)
express$log.TPM[1:floor(0.1*nrow(express))] -> Top10Pgene
sequencedata <- mclapply(mc.cores=20,1:nrow(nameseq),function(x){
  seq <- substring(as.character(nameseq$myseq[x]),seq(1,(nchar(nameseq$myseq[x])-2),by=3),seq(3,nchar(nameseq$myseq[x]),by=3))
  a <- as.data.frame(table(seq));
  names(a) <- c("codon","fre");
  a$codon <- as.vector(a$codon)
  a$gene <- nameseq$bname[x]
  a$amino <- codon$aa[match(a$codon,codon$codon)]
  a;
}) %>% rbind.fill() %>% dplyr::filter(gene %in% Top10Pgene)
sequencedata$exp <- express$average[match(sequencedata$gene,express$log.TPM)]

xijofD <- sequencedata %>%
  group_by(codon,amino) %>% dplyr::summarise(freq=sum(fre*exp)) %>% as.data.frame() %>% 
  dplyr::filter(amino != "*" & amino != "W" & amino != "M") %>%
  group_by(amino) %>% dplyr::mutate(nn=sum(freq)) %>% dplyr::mutate(xij=freq/nn)

save(xijofD,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")


####calculate cai,tai,d

##load genome seq
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")

load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")

DPcaitaiG <- mclapply(mc.cores = 10,1:nrow(nameseq),function(x){
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
    group_by(amino) %>% dplyr::summarize(Di=(sum((yij-xij)^2)^0.5)) %>% dplyr::summarize(DP=round(prod(Di)^(1/length(Di)),3)) %>% cbind(nameseq[x,3:4])
  
}) %>% rbind.fill() %>% unique()
save(DPcaitaiG,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpofecoligene.Rdata")
##ori data
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

load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpofecoligene.Rdata");
DPcaitaiG$importance <- alldata$imp[match(toupper(DPcaitaiG$bname),alldata$B.numBers)]
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/MG1655.express.31797920.csv",stringsAsFactors = F)
express$average <- (express$control__wt_glc__1 + express$control__wt_glc__2)/2
DPcaitaiG$mRNAexp <- express$average[match(DPcaitaiG$bname,express$log.TPM)]
save(DPcaitaiG,file = "/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpofecoligene.Rdata")

####
## Prediction 1 :fig4a rho (importance ~ D) in ecoli
## overall correlation (D versus importance) is negative
## because of the effect of payoff
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpofecoligene.Rdata")

geneDP0 <- DPcaitaiG %>% arrange(mRNAexp) %>% dplyr::filter(!if_any(c("mRNAexp","importance"), is.na))
rnaThresQuant <- c(1,0.5,0.1,0.05,0.04,0.03);

##cor between D and importance, group by mRNA expression
dfcor.D.imp <- lapply(rnaThresQuant,function(x){
  thisSubset <- (geneDP0 %>% arrange(desc(mRNAexp)))[1:floor(x*nrow(geneDP0)),] %>% dplyr::filter(!if_any(c("DP","importance","mRNAexp"), is.na))
  mySize <- nrow(thisSubset);
  if(mySize < 3) {
    return(data.frame(NULL))
  } else {
    pcorObj1 <- cor.test(thisSubset$DP,thisSubset$importance,method="s",alternative="less");
    pcorObj2 <- cor.test(thisSubset$DP,thisSubset$importance,method="p",alternative="less");
    
    data.frame(
      thres = rep(paste(x*100,"%",sep = ""),2),
      n = rep(mySize,2),
      estimate = c(pcorObj1$estimate,pcorObj2$estimate),
      p.val = c(pcorObj1$p.value,pcorObj2$p.value),
      corType = factor(c("rho","r"),
                       levels=c("rho","r"),
                       labels=c(expression("Spearman's"~italic(rho) ),
                                expression("Pearson's"~italic(R)) )),
      corWith = rep(c("rna"),each=2)) %>%
      return();
  }
}) %>%
  rbind.fill() %>%
  arrange(corType);
source("~/Rfunction/style.print.R")
dfcor.D.imp$thres <- factor(dfcor.D.imp$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.D.imp %>%
  dplyr::filter(corType=="\"Spearman's\" ~ italic(rho)") %>%
  ggplot(aes(y=estimate,x=thres,fill=thres,color=p.val<0.1)) +
  geom_bar(stat="identity") +
  scale_color_manual("",values=c("TRUE"="black","FALSE"="NA"),labels=c("TRUE"="P<0.05","FALSE"="N.S.")) +
  #facet_grid(~ corType, labeller = label_parsed) +
  scale_fill_brewer("Fraction of top expressed genes\n(top expression percentile)",type="seq",
                    labels = rnaThresQuant) +
  #scale_y_continuous("Partial correlation between D and importance, controlled for mRNA expression",limits = c(-0.22,0)) +
  scale_y_continuous("rho (Importance ~ D)",limits = c(-0.2,0),breaks = c(-0.2,-0.1,0)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes\n(top expression percentile)") +
  
  style.print()+theme(legend.position = "none")


#### fig 5b rho (expression ~ D)
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/dpofecoligene.Rdata")
tmpaa <- DPcaitaiG %>% dplyr::filter(!is.na(DP) & !is.na(mRNAexp))

rnaThresQuant <- c(1,0.5,0.1,0.05,0.04,0.03);
geneDP0 <- DPcaitaiG %>% arrange(mRNAexp) %>% dplyr::filter(!if_any(c("mRNAexp"), is.na))
##cor between D and importance, group by mRNA expression
dfcor.exp.D <- lapply(rnaThresQuant,function(x){
  thisSubset <- (geneDP0 %>% arrange(desc(mRNAexp)))[1:floor(x*nrow(geneDP0)),] %>% dplyr::filter(!if_any(c("DP","mRNAexp"), is.na))
  mySize <- nrow(thisSubset);
  if(mySize < 3) {
    return(data.frame(NULL))
  } else {
    pcorObj1 <- cor.test(thisSubset$DP,thisSubset$mRNAexp,method="s",alternative="greater");
    pcorObj2 <- cor.test(thisSubset$DP,thisSubset$mRNAexp,method="p",alternative="greater");
    
    data.frame(
      thres = rep(paste(x*100,"%",sep = ""),2),
      n = rep(mySize,2),
      estimate = c(pcorObj1$estimate,pcorObj2$estimate),
      p.val = c(pcorObj1$p.value,pcorObj2$p.value),
      corType = factor(c("rho","r"),
                       levels=c("rho","r"),
                       labels=c(expression("Spearman's"~italic(rho) ),
                                expression("Pearson's"~italic(R)) )),
      corWith = rep(c("rna"),each=2)) %>%
      return();
  }
}) %>%
  rbind.fill() %>%
  arrange(corType);
source("~/Rfunction/style.print.R")
dfcor.exp.D$thres <- factor(dfcor.exp.D$thres,levels = paste(rnaThresQuant*100,"%",sep = ""))
dfcor.exp.D %>%
  dplyr::filter(corType=="\"Spearman's\" ~ italic(rho)") %>%
  ggplot(aes(y=estimate,x=thres,fill=thres,color=p.val<0.1)) +
  geom_bar(stat="identity") +
  scale_color_manual("",values=c("TRUE"="black","FALSE"="NA"),labels=c("TRUE"="P<0.05","FALSE"="N.S.")) +
  scale_fill_brewer("Fraction of top expressed genes\n(top expression percentile)",type="seq",
                    labels = rnaThresQuant) +
  scale_y_continuous("rho (Expression ~ D)",limits = c(-0.02,0.35),breaks = c(0,0.1,0.2,0.3)) +
  theme_classic() +
  scale_x_discrete("Fraction of top expressed genes\n(top expression percentile)") +
  style.print()+theme(legend.position = "none")





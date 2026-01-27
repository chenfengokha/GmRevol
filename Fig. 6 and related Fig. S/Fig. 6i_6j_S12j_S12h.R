library(plyr)
library(dplyr)
library(ggplot2)
library(parallel)
library(Biostrings)
library(stringr)
library(tidyr)

functionfodpcaitai <- function(x){
  a <- substring(x,seq(1,(nchar(x)-2),by=3),seq(3,nchar(x),by=3)) %>% table() %>% as.data.frame()
  names(a)[1] <- c("codon")
  a$codon <- as.vector(a$codon)
  a$amino <- xijofD$amino[match(a$codon,xijofD$codon)]
  yima <- sum((a %>% dplyr::filter(codon %in% c("TGA","TAG","TAA")))$Freq)
  if(yima==1){
    b <- a %>% dplyr::filter(amino != "NA" & amino != "W" & amino != "M") %>%
      group_by(amino) %>% dplyr::mutate(nn=sum(Freq)) %>% dplyr::mutate(yij=Freq/nn)
    
    xijdata <- xijofD
    xijdata$yij <- b$yij[match(xijdata$codon,b$codon)]
    xijdata$yij[which(is.na(xijdata$yij))] <- 0
    resfinal <- xijdata %>% group_by(amino) %>% 
      dplyr::summarize(Diexp=(sum((yij-xij)^2))^0.5,Ditai=(sum((yij-xijoftai)^2))^0.5) %>%
      dplyr::summarize(Dpexp=prod(Diexp)^(1/length(Diexp)),Dptai=prod(Ditai)^(1/length(Ditai))) %>%
      as.data.frame()
    
  } else {
    resfinal <- data.frame(stringsAsFactors = F,Dpexp=9999999,Dptai=9999999)}
  return(resfinal)
}

#load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastcodonsupply.Rdata")
load("/mnt/data5/disk/chenfeng/NC2025review2nd/tRNAcounts/testfile/01.xijofD.SC.Rdata")

xijofD <- xijofD[,c(1:3,5)]
names(xijofD)[4] <- "xijoftai"
load("~/project/HGT/exp/zhangNature/21geneseqyeast.Rdata")

naturepaperfile <- system("ls ~/project/HGT/exp/zhangNature/*txt",intern=T)
dpcaitaiofnature <- mclapply(mc.cores = 5,1:21,function(x){
  mydf <- (read.table(naturepaperfile[x],sep = "",header = T,stringsAsFactors = F))[,1:6] %>% 
    separate(Genotype,c("site","ref","mut"),convert = T) %>% arrange(site)
  gene <- strsplit(strsplit(naturepaperfile[x],"/home/chenfeng/project/HGT/exp/zhangNature/")[[1]][2],"_all_fitness_values.txt")[[1]][1]
  myseq <- (genename %>% dplyr::filter(genename2 == gene))$myseq
  refvalue <- functionfodpcaitai(myseq)
  
  tmprefseq <- mydf[,1:2] %>% unique() 
  tmpdata <- data.frame(stringsAsFactors = F,site=1:nchar(myseq))
  tmpdata$ref <- tmprefseq$ref[match(tmpdata$site,tmprefseq$site)]
  NAsite <- which(is.na(tmpdata$ref))
  
  origseq <- myseq
  tmpsitdata <- data.frame(stringsAsFactors = F,
                           mydelta=NAsite-c(0,NAsite[-length(NAsite)]),
                           end=NAsite,sta=c(0,NAsite[-length(NAsite)])) %>%
    dplyr::filter(mydelta==max(mydelta))
  findseq <- paste((tmpdata %>% dplyr::filter(site>=(tmpsitdata$sta+1) & site<=(tmpsitdata$end-1)))$ref,collapse = "")
  match_positions <- gregexpr(findseq, origseq)[[1]][1]
  
  if(length(match_positions)>1){
    ressite <- data.frame(stringsAsFactors = F,siteonrefseq=9999999,siteonmutseq=9999999,gene)
  } else {
    ressite <- data.frame(stringsAsFactors = F,siteonrefseq=match_positions,siteonmutseq=tmpsitdata$sta+1,gene)
  }
  
  
  mclapply(mc.cores=1,1:nrow(mydf), function(y){
    newtmpmut <- mydf[y,]
    sitonrefseq <- newtmpmut$site - ressite$siteonmutseq + ressite$siteonrefseq
    
    if(substr(myseq,sitonrefseq,sitonrefseq)==newtmpmut$ref){
      newseq <- myseq
      str_sub(newseq,sitonrefseq,sitonrefseq) <- newtmpmut$mut
      mutvalue <- functionfodpcaitai(newseq) %>% cbind(newtmpmut)
    } else{
      mutvalue <- data.frame(stringsAsFactors = F,D=9999999) %>% cbind(newtmpmut)
    }
    mutvalue %>% cbind(data.frame(stringsAsFactors = F,refDexp=refvalue$Dpexp,refDtai=refvalue$Dptai,gene))
  }) %>% rbind.fill()
  
  
}) %>% rbind.fill()
dpcaitaiofnature <- dpcaitaiofnature %>% dplyr::filter(Mutation_type=="Synonymous_mutation")
dpcaitaiofnature$fitness <- (dpcaitaiofnature$Fitness_from_YPD_replicate_1+dpcaitaiofnature$Fitness_from_YPD_replicate_2+dpcaitaiofnature$Fitness_from_YPD_replicate_3+dpcaitaiofnature$Fitness_from_YPD_replicate_4)/4
save(dpcaitaiofnature,file = "/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/dpcaitaigof21gene.Rdata")

#finaldata D increase or decrease group
load("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/dpcaitaigof21gene.Rdata")
express <- (read.csv("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/21geneinfor.csv",stringsAsFactors = F))[,1:3] %>%
  arrange((Expression_level_RPKM))
dpcaitaiofnature$explevel <- express$Expression_level_RPKM[match(dpcaitaiofnature$gene,express$Gene.name)]

dpcaitaiofnature$delterType <- "F"
dpcaitaiofnature$delterType[which((dpcaitaiofnature$Fitness_from_YPD_replicate_1<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_2<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_3<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_4<1))] <- "T"

names(dpcaitaiofnature)[c(1,11)] <- c("D","refD")

##Fig.S11
Dchange <- dpcaitaiofnature %>% dplyr::mutate(deltaD=D-refD)
source("~/Rfunction/style.print.R")
Dchange$dtype <- "D-decreasing mutations"
Dchange$dtype[which(Dchange$deltaD>0)] <- "D-increasing mutations"
Dchange$gene <- factor(Dchange$gene,levels = express$Gene.name)
Dchange %>%
  ggplot(aes(gene,deltaD)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.colour = NA) +
  stat_summary(fun = 'mean', geom = 'point', shape = 18, colour = 'black')+
  facet_grid(~dtype)+
  labs(y="Changes in D",x="Genes (in ascending order of expression)")+
  style.print() +
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),legend.position =" none ")
write.csv(Dchange,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S11.csv")
###fig6i and 12g
names(dpcaitaiofnature)[c(1,2,11,12)] <- c("a","D","b","refD")
delteriousDinorde <- mclapply(mc.cores = 4,1:length(express$Gene.name),function(x){
  geneset <- express$Gene.name[x]
  mydfa <- dpcaitaiofnature %>% dplyr::filter(gene == geneset)
  
  increasD <- mydfa %>% dplyr::filter(D > refD) 
  decreasD <- mydfa %>% dplyr::filter(D < refD)
  inD1 <- wilcox.test(increasD$fitness,1,alternative="less")$p.value
  deD1 <- wilcox.test(decreasD$fitness,1,alternative="less")$p.value
  indeD <- wilcox.test(decreasD$fitness,increasD$fitness,alternative="less")$p.value
  inDdelper <- nrow(increasD %>% dplyr::filter(delterType=="T"))/nrow(increasD)
  inDdelperP <- binom.test(nrow(increasD %>% dplyr::filter(delterType=="T")),nrow(increasD),alternative = "greater",p=0.5)$p.value # 0.7751553 average deleterious mutation percentage of 21 gene
  deDdelper <- nrow(decreasD %>% dplyr::filter(delterType=="T"))/nrow(decreasD)
  deDdelperP <- binom.test(nrow(decreasD %>% dplyr::filter(delterType=="T")),nrow(decreasD),alternative = "greater",p=0.5)$p.value # 0.8661640
  data.frame(stringsAsFactors = F,percentage_deletious=c(inDdelper,deDdelper),lessthan1fit=c(inD1,deD1),morethan.5fra=c(inDdelperP,deDdelperP),type1=c("D−decreasing mutations","D−increasing mutations"),type2=c("D","D"),exp=rep(geneset,2)) %>%
    cbind(data.frame(stringsAsFactors = F,geneset,meanFit=c(mean(increasD$fitness),mean(decreasD$fitness)),sefit=c(sd(increasD$fitness)/(nrow(increasD)^0.5),sd(decreasD$fitness)/(nrow(decreasD)^0.5)),DelessthanIN=rep(indeD,2)))
  
}) %>% rbind.fill() 
source("~/Rfunction/style.print.R")
#delteriousDinorde$exp <- factor(delteriousDinorde$exp,levels = c(paste(rnaThresQuant*100,"%",sep = "")))
delteriousDinorde$morethan.5fratype <- "Nonsignificant"
delteriousDinorde$morethan.5fratype[which(delteriousDinorde$morethan.5fra>0.05 & delteriousDinorde$morethan.5fra<0.1)] <- "P < 0.1"
delteriousDinorde$morethan.5fratype[which(delteriousDinorde$morethan.5fra>0.01 & delteriousDinorde$morethan.5fra<0.05)] <- "P < 0.05"
delteriousDinorde$morethan.5fratype[which(delteriousDinorde$morethan.5fra<0.01)] <- "P < 0.01"

delteriousDinorde$morethan.5fratype <- factor(delteriousDinorde$morethan.5fratype,levels = unique(delteriousDinorde$morethan.5fratype))
color <- data.frame(stringsAsFactors = F,colortype=c("Nonsignificant","P < 0.1","P < 0.05","P < 0.01"),color=c("grey","#6BAED6","#3182BD","#08519C"))
tmpcolor <- data.frame(colortype=levels(delteriousDinorde$morethan.5fratype),stringsAsFactors = F)
tmpcolor$corlor <- color$color[match(tmpcolor$colortype,color$colortype)]

delteriousDinorde$exp <- factor(delteriousDinorde$exp,levels = express$Gene.name)

delteriousDinorde %>%
  dplyr::filter(type2=="D") %>%
  ggplot(aes(y=percentage_deletious,x=exp,fill=morethan.5fratype)) +
  geom_bar(stat="identity") +
  facet_grid(. ~ type1) +
  scale_fill_manual(values=c(tmpcolor$corlor))+
  scale_y_continuous("Fraction of deleterious mutations",limits = c(0,1.1),breaks = c(0,0.5,1)) +
  labs(x="Genes (in ascending order of expression)") +
  style.print()+theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),legend.position =" none ")
write.csv(delteriousDinorde,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S12g.csv")

##fig.6j and 12h
dpcaitaiofnature$typeD <- "D−increasing mutations"
dpcaitaiofnature$typeD[which(dpcaitaiofnature$D < dpcaitaiofnature$refD)] <- "D−decreasing mutations"

express <- (read.csv("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/21geneinfor.csv",stringsAsFactors = F))[,1:3] %>%
  arrange((Expression_level_RPKM))
#dpcaitaiofnature$explevel <- express$Expression_level_RPKM[match(dpcaitaiofnature$gene,express$Gene.name)]
express$rank <- 1:21


dpcaitaiofnature$exp <- express$rank[match(as.vector(dpcaitaiofnature$gene),express$Gene.name)]
#dpcaitaiofnature$exp <- as.vector(delteriousDinorde$exp)[match(as.vector(dpcaitaiofnature$gene),delteriousDinorde$geneset)]
dpcaitaiofnature$meanfit <- delteriousDinorde$meanFit[match(as.vector(dpcaitaiofnature$gene),delteriousDinorde$geneset)]
dpcaitaiofnature$sefit <- delteriousDinorde$sefit[match(as.vector(dpcaitaiofnature$gene),delteriousDinorde$geneset)]
dpcaitaiofnature$gene <- factor(dpcaitaiofnature$gene,levels = express$Gene.name)
dpcaitaiofnature %>%
  ggplot()+
  geom_point(aes(x=gene,y=fitness),color="grey",size=0.5)+
  geom_smooth(aes(x=exp,y=fitness),method = "lm",se=F)+
  geom_hline(yintercept = 1,color="red",linetype = "dashed")+
  geom_point(aes(x=exp,y=meanfit))+
  geom_errorbar(aes(x=exp,y=meanfit,ymax=meanfit+sefit,ymin=meanfit-sefit),color="black",width=0.4)+
  facet_grid(. ~ typeD)+
  labs(x="Genes (in ascending order of expression)",y="Fitness") +
  style.print()+
  theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1))
write.csv(dpcaitaiofnature,file = "/home/chenfeng/project/chenfengdata5/nc2025finalcheck/Fig.S12h.csv")
dpcaitaiofnature %>% group_by(typeD) %>%
  dplyr::summarize(rho=cor.test(meanfit,explevel,method = "s")$estimate,
                   p=cor.test(meanfit,explevel,method = "s")$p.value)




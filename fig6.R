
###fig 6c-f, take E. coli as an example
#1)calculate DP for all gene with one mutation (simulation)
######load refinf
ecolicds <- readDNAStringSet("~/project/HGT/exp/paper/EcoliCDS.fasta")
#genename and cds transform
nameseq <- mclapply(mc.cores=20,1:length(ecolicds),function(x){
  myname1 <- strsplit(names(ecolicds)[x]," ")[[1]][2]
  myname <- substr(myname1,7,(nchar(myname1)-1))
  bname1 <- strsplit(names(ecolicds)[x]," ")[[1]][3]
  bname <- substr(bname1,12,16)
  myseq <- as.character(ecolicds[[x]])
  data.frame(stringsAsFactors = F,myname,myseq,bname,len=nchar(myseq))
  
}) %>% rbind.fill() 
save(nameseq,file = "/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
write.csv(nameseq,file = "/mnt/data/home/chenfeng/project/HGT/exp/model/test11.csv")

#####load xij
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
tmp <- read.csv("~/project/HGT/exp/ecolicodonusage.csv",stringsAsFactors = F)
xij <- tmp %>% dplyr::filter(!(aa %in% c(NA,"M","W"))) %>% group_by(aa) %>% dplyr::mutate(n=sum(CFHEGa)) %>% group_by(aa,Codon) %>% dplyr::summarize(xij=CFHEGa/n)

dpofref <- mclapply(mc.cores=6,4001:length(ecolicds),function(x){
  myname1 <- (strsplit(names(ecolicds)[x]," ")[[1]][2])
  myname <- substr(myname1,7,(nchar(myname1)-1))
  myseq <- as.character(ecolicds[[x]])
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
      
      xijdata <- xij
      xijdata$yij <- b$yij[match(xijdata$Codon,b$codon)]
      xijdata$yij[which(is.na(xijdata$yij))] <- 0
      res <- xijdata %>% group_by(aa) %>% 
        dplyr::summarize(Di=(sum((yij-xij)^2))^0.5) %>%
        dplyr::summarize(dp=prod(Di)^(1/length(Di))) %>% 
        cbind(data.frame(stringsAsFactors = F,datatype))
      
    } else {
      res <- data.frame(stringsAsFactors = F,dp=99,datatype)
    }
    
    res
  }) %>% rbind.fill() %>% dplyr::filter(dp!=99)
  
  if(length(unique(alldp$datatype))>1){
    finres <- data.frame(stringsAsFactors = F,mutdp=alldp$dp[which(alldp$datatype=="mut")],gene=myname,wtdp=alldp$dp[which(alldp$datatype=="wt")])
    
  } else {
    finres <- data.frame(stringsAsFactors = F,mutdp=99,gene=myname,wtdp=99)
  }
  
  finres 
  
}) %>% rbind.fill() %>% dplyr::filter(wtdp!=99)


save(dpofref,file = "/mnt/data/home/chenfeng/project/HGT/exp/model/dpofrefgeneEcoli.randommutation.Rdata")

#2)calculate DP for MA genes

#filter mutant gene information
MAmutant1 <- read.csv("~/project/HGT/exp/model/MG1655MA.csv",stringsAsFactors = F) %>% dplyr::filter(type %in% c("WT","MMR"))

#####load xij
codon <- read.table("~/codonpaper/expriment/fcs/codon.txt",header = TRUE,stringsAsFactors = F)
tmp <- read.csv("~/project/HGT/exp/ecolicodonusage.csv",stringsAsFactors = F)
xij <- tmp %>% dplyr::filter(!(aa %in% c(NA,"M","W"))) %>% group_by(aa) %>% dplyr::mutate(n=sum(CFHEGa)) %>% group_by(aa,Codon) %>% dplyr::summarize(xij=CFHEGa/n)

##load name and seq of ref
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")

dpofmut <- mclapply(mc.cores=40,1:nrow(MAmutant1),function(x){
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
  xijdata <- xij
  xijdata$yij <- yijMA$yij[match(xijdata$Codon,yijMA$codon)]
  xijdata$yij[which(is.na(xijdata$yij))] <- 0
  xijdata %>% group_by(aa) %>% 
    dplyr::summarize(Di=(sum((yij-xij)^2))^0.5) %>%
    dplyr::summarize(dp=prod(Di)^(1/length(Di))) %>% 
    cbind(data.frame(stringsAsFactors = F,gene=atmp$Gene,bname=atmp$bname,pos=atmp$Position))
  
}) %>% rbind.fill()

save(dpofmut,file = "/mnt/data/home/chenfeng/project/HGT/exp/model/dpofmutgene.MAline.MG1655.Rdata")


# #2)

load("/mnt/data/home/chenfeng/project/HGT/exp/model/dpofrefgeneEcoli.randommutation.Rdata")
load("/mnt/data/home/chenfeng/project/HGT/exp/model/dpofmutgene.MAline.MG1655.Rdata")
names(dpofmut)[1] <- "mutdp"
load("/mnt/data/home/chenfeng/project/HGT/exp/model/Ecoli.nameandseq.Rdata")
dpofref$bname <- nameseq$bname[match(dpofref$gene,nameseq$myname)]
dpofmut$wtdp <- dpofref$wtdp[match(dpofmut$bname,dpofref$bname)]
express <- read.csv("/mnt/data/home/chenfeng/project/HGT/exp/model/GSE148719_Ecoli_expression.csv",stringsAsFactors = F) %>%
  arrange(desc(average)) %>% dplyr::filter(Genes %in% dpofref$bname & Genes %in% dpofmut$bname)


i=0.1
##111)top 10% highly expressed genes
experimentG <- express$Genes[1:floor(i*nrow(express))]

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
# controlG <- 
tmpgene <- data.frame(stringsAsFactors = F,type="Experiment",gene=setdiff(experimentG,controlG$bname)) %>% 
  rbind(data.frame(stringsAsFactors = F,type="Control",gene=setdiff(controlG$bname,experimentG)))
#lapply(1:2, function(k){
#k was set for top_10%_highly_expressed_genes (k=1) and low_10%_D_genes (k=2)
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
  
  subdata %>% rbind(cbind(amutdata[,c(1,2,3,5)],data.frame(stringsAsFactors = F,x=10001,type="mut")))
  
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
  #scale_x_continuous(breaks = c(5,15,25))+
  #labs(x="Percentage of D-increasing mutations among top 10% highly expressed genes",y="Frequency")+
  labs(x="Percentage of D−increasing mutations\namong genes with 10% lowest D values",y="Frequency")+
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





##fig 6i and j
functionfodpcaitai <- function(x){
  a <- substring(x,seq(1,(nchar(x)-2),by=3),seq(3,nchar(x),by=3)) %>% table() %>% as.data.frame()
  names(a)[1] <- c("codon")
  a$codon <- as.vector(a$codon)
  a$amino <- xijofD$amino[match(a$codon,xijofD$codon)]
  yima <- sum((a %>% dplyr::filter(codon %in% c("TGA","TAG","TAA")))$Freq)
  if(yima==1){
    b <- a %>% dplyr::filter(amino != "NA" & amino != "W") %>%
      group_by(amino) %>% dplyr::mutate(nn=sum(Freq)) %>% dplyr::mutate(yij=Freq/nn)
    
    xijdata <- xijofD
    xijdata$yij <- b$yij[match(xijdata$codon,b$codon)]
    xijdata$yij[which(is.na(xijdata$yij))] <- 0
    DP <- xijdata %>% group_by(amino) %>% 
      dplyr::summarize(Di=(sum((yij-xij)^2))^0.5) %>%
      dplyr::summarize(dp=prod(Di)^(1/length(Di)))
    
    resfinal <- data.frame(stringsAsFactors = F,D=DP$dp)
  } else {
    resfinal <- data.frame(stringsAsFactors = F,D=9999999)}
  return(resfinal)
}

load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/yeastcodonsupply.Rdata")
load("~/project/HGT/exp/zhangNature/21geneseqyeast.Rdata")

naturepaperfile <- system("ls ~/project/HGT/exp/zhangNature/*txt",intern=T)
dpcaitaiofnature <- mclapply(mc.cores = 5,1:21,function(x){
  mydf <- (read.table(naturepaperfile[x],sep = "",header = T,stringsAsFactors = F) %>% dplyr::filter(Mutation_type != "Nonsense_mutation"))[,1:6] %>% 
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
  
  
  mclapply(mc.cores=7,1:nrow(mydf), function(y){
    newtmpmut <- mydf[y,]
    sitonrefseq <- newtmpmut$site - ressite$siteonmutseq + ressite$siteonrefseq
    
    if(substr(myseq,sitonrefseq,sitonrefseq)==newtmpmut$ref){
      newseq <- myseq
      str_sub(newseq,sitonrefseq,sitonrefseq) <- newtmpmut$mut
      mutvalue <- functionfodpcaitai(newseq) %>% cbind(newtmpmut)
    } else{
      mutvalue <- data.frame(stringsAsFactors = F,D=9999999) %>% cbind(newtmpmut)
    }
    mutvalue %>% cbind(data.frame(stringsAsFactors = F,refD=refvalue$D,gene))
  }) %>% rbind.fill()
  
  
}) %>% rbind.fill()
dpcaitaiofnature <- dpcaitaiofnature %>% dplyr::filter(Mutation_type=="Synonymous_mutation")
dpcaitaiofnature$fitness <- (dpcaitaiofnature$Fitness_from_YPD_replicate_1+dpcaitaiofnature$Fitness_from_YPD_replicate_2+dpcaitaiofnature$Fitness_from_YPD_replicate_3+dpcaitaiofnature$Fitness_from_YPD_replicate_4)/4
save(dpcaitaiofnature,file = "/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/dpcaitaigof21gene.Rdata")

###delterious in D increase or decrease group
load("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/dpcaitaigof21gene.Rdata")
express <- (read.csv("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/21geneinfor.csv",stringsAsFactors = F))[,1:3] %>%
  arrange((Expression_level_RPKM))
dpcaitaiofnature$explevel <- express$Expression_level_RPKM[match(dpcaitaiofnature$gene,express$Gene.name)]

dpcaitaiofnature$delterType <- "F"
dpcaitaiofnature$delterType[which((dpcaitaiofnature$Fitness_from_YPD_replicate_1<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_2<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_3<1) & (dpcaitaiofnature$Fitness_from_YPD_replicate_4<1))] <- "T"

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
  data.frame(stringsAsFactors = F,percentage_deletious=c(inDdelper,deDdelper),lessthan1fit=c(inD1,deD1),morethan.5fra=c(inDdelperP,deDdelperP),type1=c("increas","decrease"),type2=c("D","D"),exp=rep(geneset,2)) %>%
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
  geom_hline(yintercept = c(0.7751553,0.8661640))+
  #scale_color_manual("",values=c("TRUE"="black","FALSE"="NA"),labels=c("TRUE"="P<0.1","FALSE"="N.S.")) +
  facet_grid(. ~ type1, labeller = label_parsed) +
  scale_fill_manual(values=c(tmpcolor$corlor))+
  #scale_fill_brewer("Fraction of top expressed genes (top expression percentile)",type="seq") +
  #scale_y_continuous("Partial correlation between D and importance, controlled for mRNA expression",limits = c(-0.22,0)) +
  scale_y_continuous("Fraction of deleterious mutations",limits = c(0,1.1),breaks = c(0,0.5,1)) +
  labs(x="Genes (in ascending order of expression)") +
  style.print()+theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1),legend.position =" none ")



dpcaitaiofnature$typeD <- "D−increasing mutations"
dpcaitaiofnature$typeD[which(dpcaitaiofnature$D < dpcaitaiofnature$refD)] <- "D−decreasing mutations"

express <- (read.csv("/mnt/data/home/chenfeng/project/HGT/exp/zhangNature/21geneinfor.csv",stringsAsFactors = F))[,1:3] %>%
  arrange((Expression_level_RPKM))
express$rank <- 1:21


dpcaitaiofnature$exp <- express$rank[match(as.vector(dpcaitaiofnature$gene),express$Gene.name)]
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


dpcaitaiofnature %>% group_by(typeD) %>%
  dplyr::summarize(rho=cor.test(meanfit,explevel,method = "s")$estimate,
                   p=cor.test(meanfit,explevel,method = "s")$p.value)

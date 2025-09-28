###5.4 TA of each gene
##5.4.1 TA fold change
ribofile <-  system("ls /home/chenfeng/project/HGT/exp/riboseq/4.ribo-seq/5.riboparser/04.density/*_rpf.txt",intern=T)
sampleid <- data.frame(stringsAsFactors = F,id=c("4_2","14_3","38_1","44_2","WT"),dp=c(0.278,0.127,0.778,0.878,99))

myTA <- mclapply(mc.cores = 3,1:length(ribofile), function(x){
  myribo <- read.csv(ribofile[x],stringsAsFactors = F,sep = "\t")
  mysample <- strsplit(strsplit(ribofile[x],".fastq.gz_rpf.txt")[[1]],"/")[[1]][11]
  mydp <- sampleid$dp[which(sampleid$id==mysample)]
  names(myribo)[7:9] <-c("f0","f1","f2")
  myribo$f0nor <- myribo$f0 
  mydata <- myribo %>% dplyr::filter(region == "cds") %>%
    group_by(name) %>% 
    dplyr::filter(max(from_tis) > 70 & from_tis > 20 & from_tis < (max(from_tis)-20)) %>%
    dplyr::summarize(TA=mean(f0nor,trim=0.2,na.rm=T)) %>%
    as.data.frame() %>% arrange(desc(TA))
  mydata$rank <- 1:nrow(mydata)
  fc <- ((mydata %>% dplyr::filter(substr(name,1,3)=="GmR"))$TA) / ((mydata %>% dplyr::filter(substr(name,1,3)!="GmR"))$TA %>% mean())
  if(length(fc)==0){
    res <- data.frame(fc=0,mysample,mydp,stringsAsFactors = F)
  } else{
    res <- data.frame(fc,mysample,mydp,stringsAsFactors = F)
  }
  res %>% cbind(mydata)
  
}) %>% rbind.fill()




###############rawdata
ribofile <-  system("ls /home/chenfeng/project/HGT/exp/riboseq/4.ribo-seq/5.riboparser/04.density/*_rpf.txt",intern=T)
sampleid <- data.frame(stringsAsFactors = F,id=c("4_2","14_3","38_1","44_2","WT"),dp=c(0.278,0.127,0.778,0.878,99))
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")

myTA <- myTA %>% group_by(mysample,mydp) %>% dplyr::mutate(mTA=mean(TA)) %>% 
  as.data.frame() %>% group_by(mysample,mydp,name,TA,mTA) %>%
  dplyr::summarize(TAratioR=TA/mTA)

TAvector <- seq(0,2,0.5)
DTrtsm <- mclapply(mc.cores = 5,1:length(ribofile), function(x){
  
  myribo <- read.csv(ribofile[x],stringsAsFactors = F,sep = "\t")
  myribo$Codon <- toupper(myribo$codon)
  
  mysample1 <- strsplit(strsplit(ribofile[x],".fastq.gz_rpf.txt")[[1]],"/")[[1]][11]
  mydp <- sampleid$dp[which(sampleid$id==mysample1)]
  names(myribo)[7:9] <-c("f0","f1","f2")
  #myribo$f0nor <- myribo$f0 + c(myribo$f1[-1],0) + c(0,myribo$f2[-length(myribo$f2)])
  myribo$f0nor <- myribo$f0 
  TAtmp <- (myTA %>% dplyr::filter(mysample==mysample1) %>% arrange(desc(TA)))
  myribo$Ta <- TAtmp$TA[match(myribo$name,TAtmp$name)]
  myribo$Taratio <- TAtmp$TAratioR[match(myribo$name,TAtmp$name)]
  
  lapply(1:length(TAvector), function(i){
    myribo %>% dplyr::filter(Ta > TAvector[i]) -> tmpdd
    tmpgene <- unique(tmpdd$name)
    myriboforcodon <- myribo %>% 
      dplyr::filter(name %in% tmpgene | substr(name,1,3) %in% c("GmR","YFP")) %>%
      #dplyr::filter(Ta > 0.09) %>% 
      group_by(name) %>% dplyr::filter(max(from_tis) > 70 & from_tis > 20 & from_tis < (max(from_tis)-20)) %>% as.data.frame() 
    
    DTofcodon <- mclapply(mc.cores = 6,1:length(unique(myriboforcodon$name)), function(y){
      tmpdf <- myriboforcodon %>% dplyr::filter(name == unique(myriboforcodon$name)[y]) %>% arrange(desc(f0nor))
      
      DTofcodonprob <- tmpdf %>% dplyr::mutate(arc=f0nor/Ta) %>%
        group_by(Codon,name) %>% 
        dplyr::summarize(DTallgene=mean(arc)) %>%
        as.data.frame() %>% dplyr::filter(!(Codon %in% c("TAA","TGA","TAG")))
      
      atmp <- tmpdf[round(0.2*nrow(tmpdf)):round(0.8*nrow(tmpdf)),]
      nlow <- atmp %>% dplyr::filter(f0nor == min(f0nor)) %>% nrow()
      lowdata0 <- (tmpdf %>% dplyr::filter(f0nor == min(atmp$f0nor)) %>% dplyr::mutate(pro=DTofcodonprob$DTallgene[match(Codon,DTofcodonprob$Codon)]) %>% arrange(desc(pro)))[1:nlow,-14]
      
      nlarge <- atmp %>% dplyr::filter(f0nor == max(f0nor)) %>% nrow()
      largedata0 <- (tmpdf %>% dplyr::filter(f0nor == max(atmp$f0nor)) %>% dplyr::mutate(pro=DTofcodonprob$DTallgene[match(Codon,DTofcodonprob$Codon)]) %>% arrange((pro)))[1:nlarge,-14]
      tmpdf %>% dplyr::filter(f0nor>min(atmp$f0nor) & f0nor<max(atmp$f0nor)) %>%
        rbind(lowdata0,largedata0) %>%
        dplyr::mutate(arc=f0nor/Ta) %>%
        group_by(Codon,name) %>%
        dplyr::summarize(DTonegene=mean(arc)) %>% as.data.frame()
    }) %>% rbind.fill()
    myDTofcodon <- DTofcodon %>%
      group_by(Codon) %>% 
      dplyr::summarize(DTallgene=mean(DTonegene)) %>%
      as.data.frame() %>% dplyr::filter(!(Codon %in% c("TAA","TGA","TAG"))) %>%
      dplyr::mutate(DTratio=DTallgene/mean(DTallgene)) %>% cbind(data.frame(type=(TAvector[i])))
    
    consump <- myribo %>% dplyr::filter(substr(name,1,3) == "GmR") %>% 
      dplyr::filter(region=="cds") %>%
      group_by(Codon) %>% dplyr::summarize(consump=sum(Taratio))
    
    myDTofcodon$consump <- consump$consump[match(myDTofcodon$Codon,consump$Codon)]
    myDTofcodon$consump[which(is.na(myDTofcodon$consump))] <- 0
    
    myDTofcodon$tai <- xijofD$tai[match(myDTofcodon$Codon,xijofD$codon)]
    
    myDTofcodon$RTSm <- myDTofcodon$consump/myDTofcodon$tai 
    
    myDTofcodon$sampleid <- mysample1
    myDTofcodon$dp <- mydp
    myDTofcodon$type2 <- paste(myDTofcodon$sampleid,myDTofcodon$dp)
    myDTofcodon
    
  }) %>% rbind.fill()
  
}) %>% rbind.fill()

##5.4.4

finalres <- mclapply(mc.cores = 11,1:length(unique(DTrtsm$type)), function(x){
  mydata <- DTrtsm %>% dplyr::filter(type == unique(DTrtsm$type)[x])
  
  mydata1 <- mydata
  
  mydata1WT <- mydata1 %>% dplyr::filter(dp==99)
  
  mydata1$WTdtratio <- mydata1WT$DTratio[match(mydata1$Codon,mydata1WT$Codon)]
  
  mydata1$dtratioR <- log2(mydata1$DTratio / mydata1$WTdtratio)
  
  pictt <- mydata1 %>% dplyr::filter(dp!=99) %>%
    dplyr::filter(dtratioR<10 & dtratioR>-10)
  
  pictt %>% as.data.frame() %>%
    group_by(dp,type) %>%
    dplyr::summarize(Rho=cor.test(consump,dtratioR,method = "s")$estimate,
                     p=cor.test(consump,dtratioR,method = "s")$p.value) %>%
    arrange(dp)
  
}) %>% rbind.fill()
save(finalres,file = "~/project/HGT/exp/riboseq/0.R/finaldata.raw.Rdata")
##################bootstrap
ribofile <-  system("ls /home/chenfeng/project/HGT/exp/riboseq/4.ribo-seq/5.riboparser/04.density/*_rpf.txt",intern=T)
sampleid <- data.frame(stringsAsFactors = F,id=c("4_2","14_3","38_1","44_2","WT"),dp=c(0.278,0.127,0.778,0.878,99))
load("/mnt/data/home/chenfeng/project/HGT/exp/fig4-5extended/ecolicodonsupply.Rdata")

myTA <- myTA %>% group_by(mysample,mydp) %>% dplyr::mutate(mTA=mean(TA)) %>% 
  as.data.frame() %>% group_by(mysample,mydp,name,TA,mTA) %>%
  dplyr::summarize(TAratioR=TA/mTA)

TAvector <- seq(0,2,0.5)
DTrtsm <- mclapply(mc.cores = 5,1:length(ribofile), function(x){
  
  myribo <- read.csv(ribofile[x],stringsAsFactors = F,sep = "\t")
  myribo$Codon <- toupper(myribo$codon)
  
  mysample1 <- strsplit(strsplit(ribofile[x],".fastq.gz_rpf.txt")[[1]],"/")[[1]][11]
  mydp <- sampleid$dp[which(sampleid$id==mysample1)]
  names(myribo)[7:9] <-c("f0","f1","f2")
  myribo$f0nor <- myribo$f0 
  TAtmp <- (myTA %>% dplyr::filter(mysample==mysample1) %>% arrange(desc(TA)))
  myribo$Ta <- TAtmp$TA[match(myribo$name,TAtmp$name)]
  myribo$Taratio <- TAtmp$TAratioR[match(myribo$name,TAtmp$name)]
  
  lapply(1:length(TAvector), function(i){
    myribo %>% dplyr::filter(Ta > TAvector[i]) -> tmpdd
    tmpgene <- unique(tmpdd$name)
    myriboforcodon <- myribo %>% 
      dplyr::filter(name %in% tmpgene | substr(name,1,3) %in% c("GmR","YFP")) %>%
      #dplyr::filter(Ta > 0.09) %>% 
      group_by(name) %>% dplyr::filter(max(from_tis) > 70 & from_tis > 20 & from_tis < (max(from_tis)-20)) %>% as.data.frame() 
    
    DTofcodon <- mclapply(mc.cores = 6,1:length(unique(myriboforcodon$name)), function(y){
      tmpdf <- myriboforcodon %>% dplyr::filter(name == unique(myriboforcodon$name)[y]) %>% arrange(desc(f0nor))
      
      DTofcodonprob <- tmpdf %>% dplyr::mutate(arc=f0nor/Ta) %>%
        group_by(Codon,name) %>% 
        dplyr::summarize(DTallgene=mean(arc)) %>%
        as.data.frame() %>% dplyr::filter(!(Codon %in% c("TAA","TGA","TAG")))
      
      atmp <- tmpdf[round(0.2*nrow(tmpdf)):round(0.8*nrow(tmpdf)),]
      nlow <- atmp %>% dplyr::filter(f0nor == min(f0nor)) %>% nrow()
      lowdata0 <- (tmpdf %>% dplyr::filter(f0nor == min(atmp$f0nor)) %>% dplyr::mutate(pro=DTofcodonprob$DTallgene[match(Codon,DTofcodonprob$Codon)]) %>% arrange(desc(pro)))[1:nlow,-14]
      
      nlarge <- atmp %>% dplyr::filter(f0nor == max(f0nor)) %>% nrow()
      largedata0 <- (tmpdf %>% dplyr::filter(f0nor == max(atmp$f0nor)) %>% dplyr::mutate(pro=DTofcodonprob$DTallgene[match(Codon,DTofcodonprob$Codon)]) %>% arrange((pro)))[1:nlarge,-14]
      tmpdf %>% dplyr::filter(f0nor>min(atmp$f0nor) & f0nor<max(atmp$f0nor)) %>%
        rbind(lowdata0,largedata0) %>%
        dplyr::mutate(arc=f0nor/Ta) %>%
        group_by(Codon,name) %>%
        dplyr::summarize(DTonegene=mean(arc)) %>% as.data.frame()
    }) %>% rbind.fill()
    myDTofcodon <- lapply(1:1000, function(j){
      genetmp <- sample(unique(DTofcodon$name),floor(0.9*length(unique(DTofcodon$name))),replace = F)
      DTofcodon %>%
        dplyr::filter(name %in% genetmp) %>%
        group_by(Codon) %>% 
        dplyr::summarize(DTallgene=mean(DTonegene)) %>%
        as.data.frame() %>% dplyr::filter(!(Codon %in% c("TAA","TGA","TAG"))) %>%
        dplyr::mutate(DTratio=DTallgene/mean(DTallgene)) %>% cbind(data.frame(rank=j,type=(TAvector[i])))
      
      
    }) %>% rbind.fill()
    
    
    consump <- myribo %>% dplyr::filter(substr(name,1,3) == "GmR") %>% 
      dplyr::filter(region=="cds") %>%
      group_by(Codon) %>% dplyr::summarize(consump=sum(Taratio))
    
    myDTofcodon$consump <- consump$consump[match(myDTofcodon$Codon,consump$Codon)]
    myDTofcodon$consump[which(is.na(myDTofcodon$consump))] <- 0
    
    myDTofcodon$tai <- xijofD$tai[match(myDTofcodon$Codon,xijofD$codon)]
    
    myDTofcodon$RTSm <- myDTofcodon$consump/myDTofcodon$tai 
    
    myDTofcodon$sampleid <- mysample1
    myDTofcodon$dp <- mydp
    myDTofcodon$type2 <- paste(myDTofcodon$sampleid,myDTofcodon$dp)
    myDTofcodon
    
    
  }) %>% rbind.fill()
  
}) %>% rbind.fill()

##5.4.4

finalres <- mclapply(mc.cores = 11,1:length(unique(DTrtsm$type)), function(x){
  mydata <- DTrtsm %>% dplyr::filter(type == unique(DTrtsm$type)[x])
  res1 <- lapply(1:1000, function(y){
    mydata1 <- mydata %>% dplyr::filter(rank==y)
    
    
    mydata1WT <- mydata1 %>% dplyr::filter(dp==99)
    
    mydata1$WTdtratio <- mydata1WT$DTratio[match(mydata1$Codon,mydata1WT$Codon)]
    
    mydata1$dtratioR <- log2(mydata1$DTratio / mydata1$WTdtratio)
    
    pictt <- mydata1 %>% dplyr::filter(dp!=99) %>%
      dplyr::filter(dtratioR<10 & dtratioR>-10)
    
    pictt %>% as.data.frame() %>%
      group_by(dp,rank,type) %>%
      dplyr::summarize(Rho=cor.test(consump,dtratioR,method = "s")$estimate,
                       p=cor.test(consump,dtratioR,method = "s")$p.value) %>%
      arrange(dp)
    
    
  }) %>% rbind.fill()
  
  
}) %>% rbind.fill()

final <- finalres %>% group_by(dp,type) %>% dplyr::summarize(m=mean(Rho),sd=sd(Rho)) %>% as.data.frame()
final$dp <- as.character(final$dp)  

save(final,file = "~/project/HGT/exp/riboseq/0.R/finaldata.Rdata")

###plot
load("~/project/HGT/exp/riboseq/0.R/finaldata.raw.Rdata")
load("~/project/HGT/exp/riboseq/0.R/finaldata.Rdata")
final$p.value <- finalres$p[match(paste(final$dp,final$type),paste(finalres$dp,finalres$type))]
final$p.type <- ">0.05"
final$p.type[which(final$p.value<0.05)] <- "<0.05"
final$p.type <- factor((final$p.type),levels = c(">0.05","<0.05"))
source("~/Rfunction/style.print.R")

final %>%
  dplyr::filter(type %in% c(0,0.5,1,1.5,2)) %>%
  ggplot(aes(x = type, y = m, linetype = dp, colour = dp)) +
  geom_line(position = position_dodge(width = 0.1)) +
  #scale_x_continuous(breaks = seq(0, 1, 0.2)) +
  geom_point(position = position_dodge(width = 0.1),aes(shape = p.type), size = 3) + 
  geom_errorbar(position = position_dodge(width = 0.1),aes(ymax = m + sd, ymin = m - sd),width=0.1) +
  scale_shape_manual(values = c(21, 16))+
  labs(x = "Cutoff of TA", y = "rho (codon decoding time relative to WT ~ codon consumption of GmR gene)") +
  style.print()
















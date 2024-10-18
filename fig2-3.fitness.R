library(plyr)
library(dplyr)
library(parallel)
library(ggplot2)
options(dplyr.summarise.inform = FALSE)
logphase <- function(x){
  tmp <- combn(nrow(x),2) %>% t() %>% as.data.frame() %>% dplyr::filter(V2-V1>1)
  tmpRP <- lapply(1:nrow(tmp), function(y){
    a <- x[tmp[y,1]:tmp[y,2],]
    t <- summary(lm(OD~time,data = a))
    data.frame(stringsAsFactors = F,DTslope=log(2)/data.frame(t$coefficients)[2,1],p=data.frame(t$coefficients)[2,4],Rsquare=t$r.squared) %>%
      cbind(a)
  }) %>% rbind.fill() %>% dplyr::filter(p<0.05) 
  mybest <- (tmpRP[,2:3] %>% unique() %>% arrange(desc(Rsquare)))[1:ceiling(nrow(tmp)*0.4),] %>% dplyr::filter(p==min(p))
  tmpRP %>% dplyr::filter(p==mybest$p & Rsquare==mybest$Rsquare)
  
}


allFiles <- system("ls ~/project/HGT/exp/fitness.epoach/GmR_*",intern=T)

alldata <- mclapply(mc.cores=10,1:length(allFiles),function(x) {
  a <- as.data.frame(read.csv(allFiles[x], skip=66,nrows=100,sep = "\t",fileEncoding="latin1",fill=T))[,c(-1,-2)]
  type <- as.numeric(strsplit(allFiles[x],"_")[[1]][2])
  a$time <- seq(from=1,to=1000,by=10)
  blank <- 0.089
  colrun <- colnames(a)[which((substr(colnames(a),1,3)== "X0."))]
  lapply(1:length(colrun), function(xx){
    data <- data.frame(OD=a[,which(colnames(a)==colrun[xx])]-blank)
    data$time <- a$time
    data$name <- colrun[xx]
    data$dp <- as.numeric(strsplit(strsplit(colrun[xx],c("_"))[[1]][1],"X")[[1]][2])
    data$seq <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][2])
    data$biorep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][3])
    data$tecrep <- as.numeric(strsplit(colrun[xx],c("_"))[[1]][4])
    data$batch <- strsplit(colrun[xx],c("_"))[[1]][5]
    if(max(data$OD)>1){
      
      data %>% dplyr::filter(OD > 0.2 & OD < 0.6) -> mydf
      mydf$OD <- log(mydf$OD)
      #log phase   
      LP <- logphase(mydf)
      
      DTvector <- log(2)/((sort(diff(LP$OD)/diff(LP$time),decreasing = T)))
      DTvector <- DTvector[which(DTvector != Inf)]
      res <- mydf[1,3:8] %>% cbind(data.frame(stringsAsFactors = F,type,DTmin1=DTvector[1],DTmin2=DTvector[2],
                                              DTmean=mean(DTvector),DTmedian=median(DTvector),DTslope=unique(LP$DTslope)))
      
    }else {
      res <- data[1,3:8] %>% cbind(data.frame(stringsAsFactors = F,type,DTmin1=99,DTmin2=99,
                                              DTmean=99,DTmedian=99,
                                              DTslope=99))
    }
    
    res
  }) %>% rbind.fill()
}) %>% rbind.fill() %>% dplyr::filter(DTmean!=99)
alldata %>% dplyr::filter(dp!=0) -> alldata
save(alldata,file = "~/project/HGT/exp/fitness.epoach/alldata.all.RPbest.Rdata")





library(tidyr)
library(ppcor)
##Fig.S9a
allFiles <- system("ls ~/project/HGT/exp/facs/CN_verify/cpn*",intern=T)

CNdata <- mclapply(mc.cores=2,1:length(allFiles),function(x) {
  type <- as.numeric(strsplit(strsplit(allFiles[x],"cpn.")[[1]][2],".csv")[[1]][1])
  a <- read.csv(allFiles[x])
  a %>% group_by(Sample) %>%
    dplyr::summarize(mCN=mean(Ratio),
                     se=sd(Ratio)/(length(Ratio)^0.5)) %>%
    cbind(type) %>% 
    separate("Sample",c("DP","seq"),"_", convert = TRUE)
  
}) %>% rbind.fill()
save(CNdata,file = "~/project/HGT/exp/review1nd/data/GmRvector.CN.Rdata")
source("~/Rfunction/style.print.R")
ggplot(data = CNdata,aes(x=DP,y=mCN))+geom_point()+geom_errorbar(aes(ymin=mCN-se,ymax=mCN+se)) +facet_wrap(~type,scales = "free")+
  geom_smooth(method = "lm",se=F)+
  labs(x="D values of GmR genes",y="Relative copy number of puc57 vector")+
  style.print()

CNdata %>% group_by(type) %>% dplyr::summarize(rho=cor.test(DP,mCN,method="s")$estimate[1],
                                               p=cor.test(DP,mCN,method="s")$p.value)

##Fig.S9b
load("/mnt/data4/disk/chenfeng/HGT/exp/facs/data/draw.expression3.Rdata")
load("~/project/HGT/exp/fig2-3.related/alldata.all.RPbest.filter.Rdata")

payoffcost <- alldata %>% 
  dplyr::mutate(dptype=paste(Dp,strain,GmCon)) %>% dplyr::filter(GmCon %in% c(0,60))

fitness <- mydf %>% group_by(dp,seq,biorep,type) %>% 
  dplyr::summarize(mfitness=mean(fitness0),sefitness=sd(fitness0)/(length(fitness0)^0.5)) %>%
  dplyr::mutate(dptype=paste(dp,seq,type))%>% dplyr::filter(type %in% c(0,60))

CNdata$dptype <- paste(CNdata$DP,CNdata$seq,CNdata$type)
pic1 <- CNdata %>% merge(payoffcost,by="dptype")
payoff.cor <- pic1 %>% group_by(type) %>% dplyr::summarize(rho=pcor.test(payoff,DP,mCN,method="s")$estimate[1],
                                             p=pcor.test(payoff,DP,mCN,method="s")$p.value) %>%
  cbind(data.frame(stringsAsFactors = F,type.name="Functional payoff"))

cost.cor <- pic1 %>% group_by(type) %>% dplyr::summarize(rho=pcor.test(cost,DP,mCN,method="s")$estimate[1],
                                             p=pcor.test(cost,DP,mCN,method="s")$p.value)%>%
  cbind(data.frame(stringsAsFactors = F,type.name="Translational cost"))

pic2 <- CNdata %>% merge(fitness,by="dptype")
#pic2 %>% ggplot(aes(DP,mfitness))+geom_point()+facet_grid(~type.x)
fit.cor <- pic2 %>% group_by(type.x) %>% dplyr::summarize(rho=pcor.test(mfitness,DP,mCN,method="s")$estimate[1],
                                               p=pcor.test(mfitness,DP,mCN,method="s")$p.value)%>%
  cbind(data.frame(stringsAsFactors = F,type.name="Fitness"))
names(fit.cor)[1] <- "type"

figs9b <- rbind(payoff.cor,cost.cor,fit.cor)
figs9b$type.name <- factor(figs9b$type.name,levels = c("Fitness","Translational cost","Functional payoff"))
figs9b$type <- factor(figs9b$type,levels = unique(figs9b$type))
figs9b$text <- NA
figs9b$text[which(figs9b$p<0.1)] <- "*"
figs9b$text[which(figs9b$p<0.05)] <- "**"
figs9b$text[which(figs9b$p<0.01)] <- "***"

figs9b$rho[which(figs9b$rho < 0)] <- figs9b$rho[which(figs9b$rho < 0)]*0.2/max(figs9b$rho[which(figs9b$rho < 0)] %>% abs())
figs9b$rho[which(figs9b$rho > 0)] <- figs9b$rho[which(figs9b$rho > 0)]*0.2/max(figs9b$rho[which(figs9b$rho > 0)] %>% abs())

source("~/Rfunction/style.print.R")
ggplot(figs9b, aes(type,type.name,fill=rho))+
  geom_tile(color="grey",alpha=2)+
  scale_fill_gradient2(low="navy",high = "firebrick3",mid="white")+
  geom_text(aes(label=text),col="black")+
  style.print()+
  theme_minimal()+
  theme(panel.border = element_rect(fill=NA,color="black", linewidth=1, linetype="solid"),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(angle=45,hjust=1, colour = 'black', size = 12),
        axis.text.y = element_text(colour = 'black', size = 12),
        plot.margin=unit(c(0.4,0.4,0.4,0.4),units=,"cm"))


#################Fig.S9c YFP-RT.PCR
library(ggplot2)
library(Matrix)
mydata <- read.csv("~/project/HGT/exp/02.yfpRTPCR/YFP_rtPCR.csv",header = T)

finalres <- lapply(1:length(unique(mydata$sample.id)), function(x){
  mya <- mydata %>% dplyr::filter(sample.id == unique(mydata$sample.id)[x])
  dxsCQ <- mya$Cq.Mean[which(mya$Gene.Name == "Genome")] %>% unique()
  
  mya %>% dplyr::filter(Gene.Name == "YFP") %>% dplyr::mutate(relativeEXP=2^(dxsCQ-Cq)) %>% group_by(D,sample.id) %>% dplyr::summarize(mEXP=mean(relativeEXP),seEXP=sd(relativeEXP)/length(relativeEXP)^0.5)
  
}) %>% rbind.fill()
source("~/Rfunction/style.print.R")

save(finalres,file = "~/project/HGT/exp/review1nd/data/YFP.RT-PCR.Rdata")
cor.test(finalres$D,finalres$mEXP,method = "s")

##################################

rank1 <- finalres %>% dplyr::filter(D < 0.33)
rank2 <- finalres %>% dplyr::filter(D > 0.33 & D < 0.66)
rank3 <- finalres %>% dplyr::filter(D > 0.66)

highlow.P <-  data.frame(stringsAsFactors = F,
                         len=c(length(rank1$D %>% unique()),length(rank2$D %>% unique()),length(rank3$D %>% unique())),
                         m=c(mean(rank1$mEXP),mean(rank2$mEXP),mean(rank3$mEXP)),
                         se=c(sd(rank1$mEXP)/(length(rank1$mEXP)^0.5),sd(rank2$mEXP)/(length(rank2$mEXP)^0.5),sd(rank3$mEXP)/(length(rank3$mEXP)^0.5)),
                         type=c("Minimum","Small","large"))


highlow.P$type <- factor(highlow.P$type,levels=c("Minimum","Small","large"))
source("~/Rfunction/style.print.R")
highlow.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(type,m,fill=type)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=m+se,ymin=m-se),position = position_dodge(width = 0.8),width = 0.3)+
  #scale_fill_manual(values = c(rep("#E6E6E6",2),rep("#B4B4B4",8),rep("#4D4D4D",3)))+
  scale_y_continuous(limits = c(0,400),breaks = c(0,200,400))+
  labs(x="D values of GmR gene",y="YFP expression level relative to dxs gene")+
  style.print()



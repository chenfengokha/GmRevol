
options(dplyr.summarise.inform = FALSE)
#############################
##1)Fig. S2

load("~/project/HGT/exp/review1nd/data.fig2a.Rdata")
source("~/Rfunction/style.print.R")
fig2a %>% ggplot(aes(dp0.1,D.trnaexp))+geom_point()+geom_vline(xintercept = c(0.33,0.66),linetype = 'dotted') + geom_hline(yintercept = c(0.36,0.53),linetype = 'dotted')+
  labs(x="D groups estimated with CUB of top 10% highly expressed genes as tRNA supply",y="D groups estimated with tRNA expression levels as tRNA supply")+
  style.print()
###################################################################################
###2 D of gmr gene. tRNA exp (GSE128812)

load("~/project/HGT/exp/fig2-3.related/alldata.all.RPbest.filter.Rdata")

highlow.P <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x]) %>% group_by(dp,seq,biorep) %>% dplyr::summarize(fitness0=mean(fitness0))
  rank1 <- a %>% dplyr::filter(dp < 0.33)
  rank2 <- a %>% dplyr::filter(dp > 0.33 & dp < 0.66)
  rank3 <- a %>% dplyr::filter(dp > 0.66)
  
  data.frame(stringsAsFactors = F,
             m=c(mean(rank1$fitness0)/mean(rank3$fitness0),mean(rank2$fitness0)/mean(rank3$fitness0),mean(rank3$fitness0)/mean(rank3$fitness0)),
             se=c(sd(rank1$fitness0/mean(rank3$fitness0))/(length(rank1$fitness0)^0.5),sd(rank2$fitness0/mean(rank3$fitness0))/(length(rank2$fitness0)^0.5),sd(rank3$fitness0/mean(rank3$fitness0))/(length(rank3$fitness0)^0.5)),
             typevalue=c("(0,1/3)","(1/3,2/3)","(2/3,1)"),
             type=c("rank1","rank2","rank3"),
             typeGmR=unique(mydf$type)[x])
  
}) %>% rbind.fill() %>% arrange(typeGmR)
highlow.P$typeGmR <- as.character(highlow.P$typeGmR)
highlow.P$typeGmR <- factor(highlow.P$typeGmR,levels=c("0","10","20","30","40","60","70","80","100","120","140","150","160","180"))
highlow.P$typevalue <- factor(highlow.P$typevalue,levels=c("(0,1/3)","(1/3,2/3)","(2/3,1)"))
source("~/Rfunction/style.print.R")
###Fig. 2e and Fig. 3f
highlow.P %>% 
  #dplyr::filter(typeGmR==120) %>%
  ggplot(aes(typeGmR,m,fill=typevalue)) +
  geom_bar(stat = "identity",position = "dodge",width = 0.8)+
  geom_errorbar(aes(ymax=m+se,ymin=m-se,color=typevalue),position = position_dodge(width = 0.8),width = 0.3)+
  scale_y_continuous(limits = c(0,1.5),breaks = c(0,0.5,1,1.5))+
  labs(x="Gentamicin concentration (ug/ml)",y="Relative fitness (relative to the large mismatch group)")+
  style.print()

### test of Fig. 2e and Fig. 3f
options(scipen = 200)
highlow.Pvalue <- lapply(1:length(unique(mydf$type)), function(x){
  a <- mydf %>% dplyr::filter(type==unique(mydf$type)[x]) %>% group_by(dp,seq,biorep) %>% dplyr::summarize(fitness0=mean(fitness0))
  rank1 <- a %>% dplyr::filter(dp < 0.33)
  rank2 <- a %>% dplyr::filter(dp > 0.33 & dp < 0.66)
  rank3 <- a %>% dplyr::filter(dp > 0.66)
  p12 <- wilcox.test(rank1$fitness0,rank2$fitness0)$p.value
  p13 <- wilcox.test(rank1$fitness0,rank3$fitness0)$p.value
  p23 <- wilcox.test(rank2$fitness0,rank3$fitness0)$p.value
  data.frame(stringsAsFactors = F,typeGmR=unique(mydf$type)[x],p12,p13,p23)
}) %>% rbind.fill() %>% arrange(typeGmR)


## file to examine auto-correlation in PCHi-C data files

##METHOD 1
## Randomly sample a bait and examine correlation between all PIRs in a given cell type

library(data.table)
library(magrittr)

DT <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab")

## get a list of baits

baits <- unique(DT$baitID)


DT[,uid:=paste(baitID,oeID,sep=':')]
setkey(DT,uid)
DT<-DT[!duplicated(uid),]
DT[,oeMid:=((oeEnd-oeStart)/2) + oeStart]

  # for a random cell type test the hypothesis that local CHiCAGO scores are more likely to be similar

cell_types <- names(DT)[16:32]
ss <- 10000000
foo<-lapply(cell_types,function(rcellt){
    byc <- lapply(1:22,function(chr){
      test <- DT[get(`rcellt`)>5 & oeChr==chr,.(oeMid,get(`rcellt`),oeChr)]
      #test <- DT[oeChr=='1',.(oeID,get(`rcellt`),oeChr)]
      idx<-do.call('cbind',split(sample.int(nrow(test),ss,replace=TRUE),1:2))
      idx<-idx[idx[,1]!=idx[,2],]
      t2<-cbind(test[idx[,1]],test[idx[,2]])
      setnames(t2,c('oe1','score1','chr1','oe2','score2','chr2'))
      ## remove sample duplicates
      t2<-unique(t2)
      t2[,c('oeDiff','scoreDiff'):=list(abs(oe1-oe2),abs(score1-score2))]
      dist <- which(t2$oeDiff<3e6)
      cbind(t2[dist,]$oeDiff,t2[dist,]$scoreDiff)
    }) %>% do.call('rbind',.)
    spe <- cor.test(byc[,1],byc[,2],method="spearman")
    -log10(spe$p.value) * sign(spe$estimate)
  }) %>% unlist

  cell_types <- gsub("\\_"," ",cell_types)
  ct <- factor(cell_types,levels=cell_types[order(foo,decreasing=TRUE)])
  plot.DT <- data.table(cell.type=ct,mlp=foo)

  library(ggplot2)
  library(cowplot)

  pp<-ggplot(plot.DT,aes(x=cell.type,y=mlp)) + geom_bar(stat="identity",fill="black") +
  theme(axis.text.x=element_text(angle = -90, hjust = 0)) + xlab("Cell Type") + ylab("-log10(P.Spearman) x sign(rho)")




## have a look at lagging

## get a dataset

lags <- c(1,2,5,10,50,100,1000)
ct<-'Total_CD4_NonActivated'
scores <- asinh(DT[baitChr==1 & get(`ct`)>5,][[ct]])

l <- length(scores)
r <-c()
par(mfrow=c(4,2))
for (j in 1:length(lags)){
       i <- lags[j]
       lagged<- scores[(1+i): l]
       laggedToo = scores[1:(l-i)]
       r[j] <- round(cor(lagged, laggedToo,method="spearman"),3)
}


cell_types <- names(DT)[16:32]
res<-lapply(cell_types,function(ct){
  #scores <- asinh(DT[baitChr==1 & get(`ct`)>5,][[ct]])
  scores <- asinh(DT[baitChr==1,][[ct]])
  l <- length(scores)
  spe <- lapply(lags,function(i){
    lagged<- scores[(1+i): l]
    laggedToo = scores[1:(l-i)]
    cte <- cor.test(lagged,laggedToo,method="spearman")
    data.table(cell.type=ct,lags=i,cor=cte$estimate,p.value=cte$p.value)
    #round(cor(lagged, laggedToo,method="spearman"),3)
  }) %>% rbindlist
  #data.table(cell.type=ct,lags=lags,cor=spe)
}) %>% rbindlist


ggplot(res,aes(x=as.factor(lags),y=cor)) + geom_bar(stat="identity",fill="white",col="black") +
facet_wrap(~cell.type)

plot.me <- res[,list(mean.cor=mean(cor),mean.mlp=mean(-log10(p.value))),by=lags]

pp<-ggplot(plot.me,aes(x=as.factor(lags),y=mean.cor)) +
geom_bar(stat="identity",fill="white",col="black") +
ylab(expression(bar(rho))) + xlab("Fragment offset (chr1)")
save_plot("~/tmp/autocorrelate.pdf",pp)

#test.chr <- split(test,test$oeChr)


test <- test.chr[['1']]
#1 randomly select a oeID


plot(abs(t2$oe1-t2$oe2),abs(t2$score1-t2$score2))


lapply(split(idx,1:2),function(x){
  test[]
})
test[idx,]


## take average across cell types

mean.DT<-DT[,cell_types,with=FALSE] %>% as.matrix %>% rowMeans %>% cbind(DT[,.(baitID,oeID)],meanScore=.) %>% unique
## next remove interactions with below median number of interactions
sum.DT <- mean.DT[,list(n=.N),by=baitID]
mean.DT <- mean.DT[baitID %in% sum.DT[n>median(n),]$baitID,]


## what is median number of interactions

getMean<-function(f.DT){
  sum.DT <- f.DT[,list(cor=abs(cor(meanScore,oeID)),n=.N),by=baitID]
  med <- median(sum.DT$n)
  mean(sum.DT[n>med,]$cor)
}

act <- getMean(mean.DT)
shuff <- copy(mean.DT)
perm <- sapply(1:100,function(i){
  message(i)
  ## for each perm we shuffle correlations - thus keep number of interactions for a bait the same
  getMean(shuff[,meanScore:=sample(meanScore)])
})

## what about if we use t statistics from cor.test without censoring by median

mean.DT<-DT[,cell_types,with=FALSE] %>% as.matrix %>% rowMeans %>% cbind(DT[,.(baitID,oeID)],meanScore=.) %>% unique
## next remove interactions with below median number of interactions
sum.DT <- mean.DT[,list(n=.N),by=baitID]
mn<-median(sum.DT$n)
mn<-2
mean.DT <- mean.DT[baitID %in% sum.DT[n>mn,]$baitID,]


## what is median number of interactions

getMean<-function(f.DT){
  sum.DT <- f.DT[,list(cor=abs(cor.test(meanScore,abs(baitID-oeID))$statistic),n=.N),by=baitID]
  mean(sum.DT$cor)
  #med <- median(sum.DT$n)
  mean(sum.DT$cor)
}

act <- getMean(mean.DT)
shuff <- copy(mean.DT)
n.perms<-100
perm <- sapply(1:n.perms,function(i){
  message(i)
  ## for each perm we shuffle correlations - thus keep number of interactions for a bait the same
  getMean(shuff[,meanScore:=sample(meanScore)])
})

P <- sum(act<perm)/n.perms

Z <- act-mean(perm)/sqrt(var(perm))

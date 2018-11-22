## analyse hierachical prioritisation

DATA_DIR <- '/home/ob219/share/cogs_bb/COGS'

fs <- list.files(path=DATA_DIR,pattern="\\_prioritised.tab",full.names=TRUE)

res.DT <- lapply(fs,fread) %>% rbindlist

## examine protein_coding only

resf <- res.DT[biotype=='protein_coding',]

## filter such that there is at least marginal overall posterior that a  gene
## is involved in a disease
D <- resf

Df <- melt(D,id.vars=c('disease','ensg','name'),measure.vars='node')
c.DT<-Df[,list(gene.count=.N),by=c('disease','value')]
c.DT[,]
test.DT <- c.DT[disease %in% c('CD','RA','SLE','T1D','UC'),]
back.DT <- c.DT[!disease %in% c('CD','RA','SLE','T1D','UC'),]

bdist <- back.DT[,list(mean=mean(log(gene.count+1)),var=var(log(gene.count+1))),by='value']

M <- merge(test.DT,bdist,by='value')
M[,Z:=(log(gene.count+1)-mean)/sqrt(var)]

## try Tukey transform

back.DT[,tukey.trans:=transformTukey(gene.count,plotit=FALSE),by='value']

bdist <- back.DT[,list(mean=mean(log(gene.count+1)),var=var(log(gene.count+1))),by='value']


setnames(M,'value','node.name')
M <- melt(M,id.vars=c('node.name','disease'),measure.vars='Z')
#M <- dcast(M,disease~node.name+variable)
#mat <- as.matrix(M)[,-1] %>% apply(.,2,as.numeric)
#rownames(mat) <- M$disease
#colnames(mat) <- colnames(mat) %>% sub("\\_Z","",.)
#library(pheatmap)
#pheatmap(mat)

## do in ggplot
M[,p.adj:=(pnorm(value,lower.tail=FALSE) * 2) %>% p.adjust]
library(cowplot)
ggplot(M[p.adj<0.05,],aes(x=node.name,y=disease,fill=value)) + geom_tile(color='black') +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + xlab("Node") + ylab("Disease") +
scale_fill_continuous("Z")


## looking at background distro

node <- 'overall'

library(rcompanion)
par(mfrow=c(2,2))
qqnorm(back.DT[value==node,]$gene.count,main="Overall node gene counts")
qqline(back.DT[value==node,]$gene.count,col='red')
qqnorm(log(back.DT[value==node,]$gene.count + 1),main="Overall node log(gene counts+1)")
qqline(log(back.DT[value==node,]$gene.count + 1),col='red')
qqnorm(sqrt(back.DT[value==node,]$gene.count),main="Overall node sqrt(gene counts)")
qqline(sqrt(back.DT[value==node,]$gene.count),col='red')
T_tuk = transformTukey(back.DT[value==node,]$gene.count,plotit=FALSE)
qqnorm(T_tuk,main="Overall node transformTukey")
qqline(T_tuk,col='red')
par(mfrow=c(1,1))
hist(back.DT[value==node,]$gene.count,main="Overall node gene counts")

bb_trait <- readRDS('/home/ob219/share/cogs_bb/bb_trait_manifest.RDS')

scount.DT <- merge(back.DT,bb_trait[,.(phe,cases)],by.y='phe',by.x='disease')

test <- scount.DT[,list(cor=cor(gene.count,cases)),by='value']

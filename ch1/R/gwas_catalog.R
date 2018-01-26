library(data.table)
library(ggplot2)
library(magrittr)

DT <- fread('https://www.ebi.ac.uk/gwas/api/search/downloads/studies')
setnames(DT,make.names(names(DT)))
## get the list of unique studies
setkey(DT,PUBMEDID)
DT.u <- unique(DT[,.(DATE,DISEASE.TRAIT)])
DT.u[,DATE:=as.Date(DATE)]
DT.u <- DT.u[order(DATE),]
DT.u[,es:=as.integer(DATE)]
sum<-DT.u[,list(studyno=.N),by=DATE]
sum[,cs:=cumsum(studyno)]
sum[,group:='group1']

library(cowplot)

p<-ggplot(sum,aes(x=DATE,y=cs)) + geom_line() +
xlab("Study Publication Date") + ylab("Cum. Study Count") +
geom_vline(xintercept=as.integer(as.Date("2007-06-04")),color='firebrick',lty=2) +
geom_vline(xintercept=as.integer(as.Date("2005-03-10")),color='dodgerblue',lty=2)


p <- p + annotate("text", x = as.Date("2007-06-04") + 70, y = 2000, label = "WTCCC (2007)",color='firebrick',angle=90) +
annotate("text", x = as.Date("2005-03-10")+ 70, y = 2000, label = "Klein et al. (2005)",color='dodgerblue',angle=90)

DT.as <-  fread('https://www.ebi.ac.uk/gwas/api/search/downloads/alternative')
setnames(DT.as,make.names(names(DT.as)))
DT.asf <- DT.as[as.numeric(P.VALUE)<5e-8 & SNP_ID_CURRENT!='' ,.(PUBMEDID,DATE,MAPPED_TRAIT,SNP_ID_CURRENT,P.VALUE,OR.or.BETA,MAPPED_TRAIT_URI)]
DT.asf[,rs:=sprintf("rs%s",DT.asf$SNP_ID_CURRENT)]

library("SNPlocs.Hsapiens.dbSNP144.GRCh37")
snps<-SNPlocs.Hsapiens.dbSNP144.GRCh37
lu<-snpsById(snps,unique(sprintf("rs%s",DT.asf$SNP_ID_CURRENT)),ifnotfound='drop')
lu<-data.table(as.data.frame(lu))
setkey(lu,'RefSNP_id')
setkey(DT.asf,'rs')
DT.asf<-DT.asf[lu]
t.gr <- with(DT.asf,GRanges(seqnames=Rle(gsub("ch","",seqnames)),ranges=IRanges(start=pos,width=1L),line=1:nrow(DT.asf)))

## next assign these to EUR haplotype block.

blocks <- fread('/home/ob219/scratch/as_basis/support_tab/all.1cM.tab')[,c('id','width'):=list(1:.N,end-start+1)]
setkey(blocks,id)
blocks.gr<-with(blocks,GRanges(seqnames=Rle(chr),ranges=IRanges(start=start,end=end),id=1:nrow(blocks)))
ol<-as.matrix(findOverlaps(t.gr,blocks.gr))
DT.asf[,ld.block:=0]
DT.asf[ol[,1],ld.block:=ol[,2]]

## next we remove those not in a blocks

DT.filt<-DT.asf[ld.block!=0,]
DT.filt<-DT.filt[,DATE:=as.Date(DATE)]
DT.filt<-DT.filt[order(DATE),]
DT.filt <- DT.filt[!duplicated(ld.block),]
setkey(DT.filt,ld.block)
DT.filt <- DT.filt[blocks][!is.na(rs),]
pmin(DT.filt$width,median(DT.filt$width))
#
sum.as<-DT.filt[,list(width=sum(pmin(width,2e6))),by=DATE]
sum.as <- sum.as[order(DATE),]
sum.as[,cs.mb:=signif(cumsum(width)/3e9,digits=3)*100]

pl<-ggplot(sum.as,aes(x=DATE,y=cs.mb)) + geom_line() +
xlab("Study Publication Date") + ylab("% of genome assoc.") +
geom_hline(yintercept=50,color='firebrick',lty=2)
ppf <- plot_grid(p, pl, labels = c("A", "B"),nrow = 2, align = "v")
save_plot("~/tmp/gwas_plot.pdf",ppf,base_height=6)


## chris asked for the rate of change of study/rate of change of region vs years

## intervals
sum.as[,elapsed:=c(0,diff(as.numeric(sum.as$DATE) - min(as.numeric(sum.as$DATE))))]
sum.as[,rate:=width/((elapsed+1)*10e6)]
pg<-ggplot(sum.as,aes(x=DATE,y=rate)) + geom_line() +
xlab("Study Publication Date") + ylab("Delta") +
geom_vline(xintercept=as.integer(as.Date("2007-06-04")),color='firebrick',lty=2) +
geom_vline(xintercept=as.integer(as.Date("2005-03-10")),color='dodgerblue',lty=2)

## something similar for sudy count

sum[,elapsed:=c(0,diff(as.numeric(sum$DATE) - min(as.numeric(sum$DATE))))]
sum[,rate:=studyno/((elapsed+1))]
ps<-ggplot(sum,aes(x=DATE,y=rate)) + geom_line() +
xlab("Study Publication Date") + ylab("Delta") +
geom_vline(xintercept=as.integer(as.Date("2007-06-04")),color='firebrick',lty=2) +
geom_vline(xintercept=as.integer(as.Date("2005-03-10")),color='dodgerblue',lty=2)

## merge
setnames(sum,'DATE','sDATE')
setnames(sum.as,'DATE','aDATE')

sum[,sDATE2:=c(sDATE[-1],max(sDATE)+2)-1]
sum[,sid:=1:.N]
sum.as[,aDATE2:=c(aDATE[-1],max(aDATE)+2)-1]
sum.as[,aid:=1:.N]
setkey(sum,sDATE,sDATE2)
setkey(sum.as,aDATE,aDATE2)
olbd <- foverlaps(sum,sum.as)

## study is not duplicated therefore sum by aid
merged <- olbd[,list(total.study=sum(studyno),DATE=as.Date(unique(aDATE)),width=unique(width)),by=aid][!is.na(DATE),]
merged[,elapsed:=c(0,diff(as.numeric(merged$DATE) - min(as.numeric(merged$DATE))))]
merged[,c('pub.rate','assoc.rate'):=list(total.study/((elapsed+1)),width/((elapsed+1)))]
merged[,overall.rate:=assoc.rate/pub.rate]

## there are three main steps what happened here ?
interest.dates <- head(merged[order(overall.rate,decreasing=TRUE),]$DATE,n=3)


ops<-ggplot(merged,aes(x=DATE,y=cumsum(overall.rate))) + geom_line() +
#xlab("Study Publication Date") + ylab("Associated recomb.block width rate/Publication Rate") +
xlab("Study Publication Date") + ylab(TeX('Cumulative $\\frac{\\Delta Basepair Associated}{\\Delta Publication}$')) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[1])),color='firebrick',lty=2) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[2])),color='dodgerblue',lty=2) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[3])),color='steelblue',lty=2) +
annotate("text", x = as.Date(interest.dates[1]) + 70, y = 4e8, label = "Height:Lettre..;Gudbjartsson..;Weedon..",color='firebrick',angle=90) +
annotate("text", x = as.Date(interest.dates[2])+ 70, y = 4e8, label = "Height:Lango-Allen et al.",color='dodgerblue',angle=90) +
annotate("text", x = as.Date(interest.dates[3]) + 70, y = 4e8, label = "Schizophrenia:Ripke et.al",color='steelblue',angle=90)

ops<-ggplot(merged,aes(x=DATE,y=overall.rate)) + geom_line() +
#xlab("Study Publication Date") + ylab("Associated recomb.block width rate/Publication Rate") +
xlab("Study Publication Date") + ylab(TeX('$\\frac{\\Delta Basepair Associated}{\\Delta Publication}$')) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[1])),color='firebrick',lty=2) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[2])),color='dodgerblue',lty=2) +
geom_vline(xintercept=as.integer(as.Date(interest.dates[3])),color='steelblue',lty=2) +
annotate("text", x = as.Date(interest.dates[1]) + 70, y = 4e7, label = "Height:Lettre..;Gudbjartsson..;Weedon..",color='firebrick',angle=90) +
annotate("text", x = as.Date(interest.dates[2])+ 70, y = 4e7, label = "Height:Lango-Allen et al.",color='dodgerblue',angle=90) +
annotate("text", x = as.Date(interest.dates[3]) + 70, y = 4e7, label = "Schizophrenia:Ripke et.al",color='steelblue',angle=90)

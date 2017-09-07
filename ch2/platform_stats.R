library(data.table)

## code to generate plaform stats for

## annotations used in Javierre et al.
DT<-fread('/scratch/ob219/DATA/JAVIERRE_GWAS/support/HindIII_baits_e75.bed')
## split V5
DT[, c('name','ensg','bt','strand'):=tstrsplit(V5,':',fixed=TRUE)]

## load in gtf file
GTF<-fread('/scratch/ob219/DATA/JAVIERRE_GWAS/support/Homo_sapiens.GRCh37.75.genes.gtf')
GTF[,ensg:=gsub('.*"(.*)"',"\\1",tstrsplit(GTF$V9,';',fixed=TRUE,fill=TRUE)[[1]])]
GTF<-GTF[GTF$V1 %in% unique(DT$V1),]
## bin ensg by biotype

gbt<-lapply(split(GTF$ensg,GTF$V2),unique)
gbtc<-sapply(gbt,length)

summ<-data.table(class=names(gbtc),overall=gbtc,key='class')
pbt<-lapply(split(DT$ensg,DT$bt),unique)
pbtc<-sapply(pbt,length)
summ.pc<-data.table(class=names(pbtc),pc=pbtc,key='class')
summ<-summ[summ.pc]
summ[,perc:=round((pc/overall)*100)]
summ<-summ[order(summ$perc,decreasing=TRUE),]

library(xtable)
xtable(summ)


## what is the distribution of non protein coding genes across the genome ?
bait.pc<-unique(subset(DT,bt=='protein_coding')$ensg)
GTF.pc.nc<-subset(GTF,V2=='protein_coding' & !ensg %in% bait.pc )

## ok get the canonical TSS for this analysis - whilst I realise
GTF.pc.nc[,coord:=V4]
GTF.pc.nc[GTF.pc.nc$V7=='-',coord:=V5]
neg<-GTF.pc.nc[GTF.pc.nc$V7=='-',]
ps<-GTF.pc.nc[GTF.pc.nc$V7=='+',]
gr<-rbind(ps[ps[, .I[which.min(V4)], by=ensg]$V1],neg[neg[, .I[which.min(V4)], by=ensg]$V1])[,.(V1,coord)]
gr<-with(gr,GRanges(seqnames=Rle(V1),range=IRanges(start=coord,end=coord)))
## get

## get max coord for each chr
library(GenomicRanges)
max.coord<-GTF.pc.nc[,list(mc=max(V5)),by=V1]
reso<-1e7
genome<-lapply(split(max.coord,max.coord$V1),function(m){
    s<-seq(1,m$mc,by=reso)
    s2<-c(head(s+reso,-1)-1,m$mc)
    GRanges(seqnames=Rle(m$V1),ranges=IRanges(start=s,end=s2))
})
genome<-unlist(GRangesList(genome))
ol<-as.matrix(findOverlaps(genome,gr))
counts<-sapply(split(ol[,2],ol[,1]),length)
genome$gc<-0
genome[as.numeric(names(counts)),]$gc<-counts
gt<-data.table(data.frame(genome))
gt$chr<-as.numeric(as.character(gt$seqnames))
gt[gt$seqnames=='X']$chr<-23
gt[gt$seqnames=='Y']$chr<-24
gt<-gt[order(gt$chr,gt$start),]

tmp<-split(gt$end,gt$chr)
tmp2<-split(gt$start,gt$chr)

cs<-c(0,head(cumsum(as.numeric(sapply(tmp,max))),-1))+1
for(i in seq_along(tmp)){
  tmp2[[i]]<-tmp2[[i]] + cs[i]
}

gt$pstart<-do.call('c',tmp2)
gt[,pend:=pstart+width]

library(ggplot2)

ggplot(gt,aes(xmin=pstart,xmax=pend,ymin=0,ymax=gc,fill=as.factor(chr%%2))) + geom_rect() + theme_bw()

ggplot(gt,aes(x=pstart,y=gc,color=as.factor(chr%%2))) + geom_path() + theme_bw()


## what is the average size of captured sequence

round(mean(width(with(DT,IRanges(start=V2,end=V3)))))

## what is the average size for all HindIII for the genome

h3<-fread('/scratch/ob219/DATA/JAVIERRE_GWAS/support/Digest_Human_HindIII.bed')
round(mean(width(with(h3,IRanges(start=V2,end=V3)))))


## protein coding promoters captured per gene
DT.pc<-DT[DT$bt=='protein_coding',]
by.frag<-split(DT.pc$ensg,DT.pc$V4)
c<-data.table(table(sapply(by.frag,function(f) length(unique(f)))))
(1-(c[1,]$N/length(by.frag))) * 100

## what about the whole shooting match



#    1     2     3     4     5     6     7
#13599  2567   347    73    13     8     1



## generate a table of coverage

## how many protein coding covered

## this is a plot that attempts to use data from delaange in order to demonstrate the utility of COGS score.
## should be run on mac

library(Gviz)
library(biomaRt)
library(data.table)
library(GenomicInteractions)

cd.DT <- fread("gunzip -c /Users/oliver/DATA/THESIS/ch3/cd_build37_40266_20161107.txt.gz")
cd.DT[,c('chr','pos','a1','a2'):=tstrsplit(MarkerName,":|_")]
#get ld block data
ld.DT <- fread("/Users/oliver/DATA/JAVIERRE_GWAS/support/1cM_regions.b37.bed")
ld.gr <- with(ld.DT,GRanges(seqnames=Rle(paste0('chr',V1)),ranges=IRanges(start=V2,end=V3),id=1:nrow(ld.DT)))


## variant of interest

vint <- cd.DT[chr=='3' & pos=='188401160',]
# gene of interest
gint <- 'BCL6'
# tissue of interest
tint <- 'Macrophages_M1'



cd.gr <- with(cd.DT,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(pos),width=1L),mlp=-log10(P.value)))

vint.gr <- GRanges(seqnames=Rle('chr3'),ranges=IRanges(start=188401160,width=1L))

## get ld block of interest

vld.gr <- subsetByOverlaps(ld.gr,vint.gr)

## next load in pm data

jav.DT <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab')
fjav.DT <- jav.DT[name==gint & get(tint)>5,]


bait.gr <- with(fjav.DT,GRanges(paste0('chr',baitChr),IRanges(start=baitStart,end=baitEnd),idx=1:nrow(fjav.DT)))
oe.gr <- with(fjav.DT,GRanges(paste0('chr',oeChr),IRanges(start=oeStart,end=oeEnd),idx=1:nrow(fjav.DT)))

## just select oe in LD block

oe.idx <- subsetByOverlaps(oe.gr,vld.gr)$idx
interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=1)

regions.gr <- reduce(subsetByOverlaps(ld.gr,GRangesList(list(bait.gr,oe.gr))),min.gapwidth=10L)
regions.gr <- GRanges(seqnames=unique(seqnames(bait.gr)),ranges=IRanges(start=min(start(c(bait.gr,oe.gr))),end=max(end(c(bait.gr,oe.gr)))))


e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

gene <- BiomartGeneRegionTrack(
    genome="hg19", name="Genes", transcriptAnnotation="symbol",
    mart=e75.genemart,
    collapseTranscripts="meta",
    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,chromosome_name='3',
    start=start(regions.gr),end=end(regions.gr),biotype='protein_coding',hgnc_symbol=list('BCL6','LPP')),
    fill="black",col="black",cex.group = 1,fontcolor.group="black",just.group="left")

tranges <- ranges(gene)
tranges <- subset(tranges,symbol  %in% c('LPP','BCL6'))

ranges(gene)<-tranges

ideo<-fread("/Users/oliver/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

## can we label mlp's by colour ?

mlp.gr <- subsetByOverlaps(cd.gr,regions.gr)
ol <- as.matrix(findOverlaps(mlp.gr,vld.gr))
mlp.gr$interest <- mlp.gr$mlp
mlp.gr[-ol[,1],]$interest <- NA
mlp.gr[ol[,1],]$mlp <- NA


tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(regions.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack(regions.gr)
tracks$gene<-gene
tracks$interactions <- InteractionTrack(name='PCHi-C\nMacrophages', interaction)
displayPars(tracks$interactions)$anchor.height <- 0.1
displayPars(tracks$interactions)$col.anchors.line <- "black"
displayPars(tracks$interactions)$col.anchors.fill <- "darkgreen"
displayPars(tracks$interactions)$col.interactions <- "darkgreen"
tracks$mut<-DataTrack(name="deLaange et al. CD\n-log10(P)",mlp.gr,groups=c('mlp','interest'),legend=FALSE,col=c('red','lightgrey'),fill=c('red','lightgrey'))
displayPars(tracks$axis)$background.title <- 'black'
displayPars(tracks$interactions)$background.title <- 'black'
displayPars(tracks$gene)$background.title <- 'black'
displayPars(tracks$mut)$background.title <- 'black'
#tracks$ld<-AnnotationTrack(name='LD',subsetByOverlaps(ld.gr,regions.gr))
#col="red",fill="red", groupAnnotation = "group",cex.group=1,fontcolor.group="red",just.group = "above")
#tracks$se <- AnnotationTrack(name="SE Hnisz",se.gr,col="orange",fill="orange")
#feature(tracks$mut) <- DT$mut.name
plotTracks(tracks)



e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

hgnc_name<-'CTLA4'

target<-BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        external_name = hgnc_name,
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T,hgnc_symbol=hgnc_name))

rchr <- gsub("^chr","",unique(seqnames(target)))
rstart<-min(start(range(target)))
rend<-max(end(range(target)))

ideo<-fread("/rds/user/ob219/hpc-work/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

ldb<-unique(subset(basis.DT,chr==rchr & between(position,rstart-1e4,rend+1e4))$ld)
raw.dat<-subset(basis.DT,ld.block %in% c(ldb-1,ldb,ldb+1))
raw.dat[,slp:=sign(log(or)) * -log10(p.value)]
melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('slp'))
slp<-dcast(melt.DT[variable=="slp",],pid+chr+position~trait+variable)
comb.gr <- with(slp,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(comb.gr)<- as.data.frame(slp[,4:ncol(slp)]) %>% DataFrame
bshrink.gr <- with(shrink.DT[pid %in% slp$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),shrink=bshrink))
wshrink.gr <- with(shrink.DT[pid %in% slp$pid,],GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(position),width=1L),wshrink=ws_ppi))
#compute Wakefields
ppi<-raw.dat[,list(pid=pid,ppi=wakefield_pp(p.value,maf,unique(n),unique(n1/n))),by=c('trait','ld.block')][,.(trait,ppi,pid,ld.block)]
setkeyv(ppi,c('pid','trait'))
setkeyv(raw.dat,c('pid','trait'))
raw.dat <- ppi[raw.dat]
melt.DT<-melt(raw.dat,id.vars=c('pid','chr','position','trait'),measure.vars=c('slp','ppi'))
ppid<-dcast(melt.DT[variable=="ppi",],pid+chr+position~trait+variable)
ppi.gr <- with(ppid,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=position,width=1L)))
mcols(ppi.gr)<- as.data.frame(ppid[,4:ncol(ppid)]) %>% DataFrame


tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(comb.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        chromosome = seqlevels(comb.gr), start = min(slp$position), end = max(slp$position),
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T))
tracks$combmlp<- DataTrack(name="-log10(P)",comb.gr,groups=gsub("_slp","",names(mcols(comb.gr))))
displayPars(tracks$combmlp)$background.title <- 'black'
tracks$ppi <- DataTrack(name="PPA",ppi.gr,groups=gsub("_ppi","",names(mcols(ppi.gr))))
tracks$shrink <- DataTrack(name="PPA Shrink",bshrink.gr,col='dodgerblue')
displayPars(tracks$shrink)$background.title <- 'dodgerblue'
tracks$shrink2 <- DataTrack(name="Weighted PPA Shrink",wshrink.gr,col='firebrick')
displayPars(tracks$shrink2)$background.title <- 'firebrick'
#plotTracks(tracks,from=min(raw.dat[ld.block %in% ldb,]$position),to=max(raw.dat[ld.block %in% ldb,]$position))

## zoom into ld block with highest posterior
tmp.DT<-shrink.DT[ld.block %in% ldb,list(tot=sum(ws_ppi)),by=ld.block]
max.ldb <- tmp.DT[which.max(tot),]$ld.block
pdf("~/tmp/ctla4_shrinkage.pdf")
plotTracks(tracks,from=min(raw.dat[ld.block %in% max.ldb,]$position),to=max(raw.dat[ld.block %in% max.ldb,]$position))
dev.off()

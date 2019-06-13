## this is a plot that attempts to use data from delaange in order to demonstrate the utility of COGS score.
## should be run on mac

library(Gviz)
library(biomaRt)
library(data.table)
library(GenomicInteractions)

cd.DT <- fread("gunzip -c /Users/oliver/DATA/THESIS/ch3/ibd_build37_59957_20161107.txt.gz")
cd.DT[,c('chr','pos','a1','a2'):=tstrsplit(MarkerName,":|_")]
#get ld block data
ld.DT <- fread("/Users/oliver/DATA/JAVIERRE_GWAS/support/1cM_regions.b37.bed")
ld.gr <- with(ld.DT,GRanges(seqnames=Rle(paste0('chr',V1)),ranges=IRanges(start=V2,end=V3),id=1:nrow(ld.DT)))


## variant of interest

vint <- cd.DT[chr=='6' & pos=='138002175',]
# gene of interest
gint <- 'TNFAIP3'
# tissue of interest
tint <- 'Total_CD4_Activated'


cd.gr <- with(cd.DT,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=as.numeric(pos),width=1L),mlp=-log10(P.value)))

vint.gr <- GRanges(seqnames=Rle('chr6'),ranges=IRanges(start=138002175,width=1L))

## get ld block of interest

vld.gr <- subsetByOverlaps(ld.gr,vint.gr)

## next load in pm data

jav.DT <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab')
fjav.DT <- jav.DT[name==gint & get(tint)>5,]


bait.gr <- with(fjav.DT,GRanges(paste0('chr',baitChr),IRanges(start=baitStart,end=baitEnd),idx=1:nrow(fjav.DT)))
oe.gr <- with(fjav.DT,GRanges(paste0('chr',oeChr),IRanges(start=oeStart,end=oeEnd),idx=1:nrow(fjav.DT)))

## just select oe in LD block

oe.idx <- subsetByOverlaps(oe.gr,vld.gr)$idx
#interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=1)

counts <- rep(1,length(oe.idx))
counts[queryHits(findOverlaps(oe.gr[oe.idx,],vint.gr))]<-2
interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=counts)



#regions.gr <- reduce(subsetByOverlaps(ld.gr,GRangesList(list(bait.gr,oe.gr))),min.gapwidth=10L)
#regions.gr <- reduce(subsetByOverlaps(ld.gr,GRangesList(list(bait.gr,oe.gr))))
regions.gr <- GRanges(seqnames=unique(seqnames(bait.gr)),ranges=IRanges(start=min(start(c(bait.gr,oe.gr))),end=max(end(c(bait.gr,oe.gr)))))


e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")

gene <- BiomartGeneRegionTrack(
    genome="hg19", name="Genes", transcriptAnnotation="symbol",
    mart=e75.genemart,
    collapseTranscripts="meta",
    stackHeight=0.2, filters=list(with_ox_refseq_mrna=T,chromosome_name='6',
    start=start(regions.gr),end=end(regions.gr),biotype='protein_coding',hgnc_symbol=list('TNFAIP3')),
    fill="black",col="black",cex.group = 1,fontcolor.group="black",just.group="left")

tranges <- ranges(gene)
tranges <- subset(tranges,symbol  %in% c('TNFAIP3'))

ranges(gene)<-tranges

ideo<-fread("/Users/oliver/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

## can we label mlp's by colour ?

mlp.gr <- subsetByOverlaps(cd.gr,regions.gr)
ol <- as.matrix(findOverlaps(mlp.gr,vint.gr))
mlp.gr$interest <- mlp.gr$mlp
mlp.gr[-ol[,1],]$interest <- NA
mlp.gr[ol[,1],]$mlp <- NA

interaction.cols <- rep('grey',length(interaction))
interaction.cols[queryHits(findOverlaps(interaction,vint.gr))]<-'red'

interaction.line.cols <- rep('grey',length(interaction))
interaction.line.cols[queryHits(findOverlaps(interaction,vint.gr))+1]<-'red'

interaction.line.cols.black <- rep('black',length(interaction))
interaction.line.cols.black[queryHits(findOverlaps(interaction,vint.gr))+1]<-'red'


tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(regions.gr),bands=ideo)
tracks$axis<-GenomeAxisTrack(regions.gr)
tracks$gene<-gene
tracks$interactions <- InteractionTrack(name='PCHi-C\nAct CD4+\n(Javierre et al.)', interaction)
displayPars(tracks$interactions)$anchor.height <- 0.1
displayPars(tracks$interactions)$col.anchors.line <- interaction.line.cols.black
displayPars(tracks$interactions)$col.anchors.fill <- interaction.line.cols
displayPars(tracks$interactions)$col.interactions <- interaction.cols
tracks$mut<-DataTrack(name="IBD-log10(P)\n(de Lange et al.)",mlp.gr,groups=c('mlp','interest'),legend=FALSE,col=c('red','lightgrey'),fill=c('red','lightgrey'))
displayPars(tracks$axis)$background.title <- 'black'
displayPars(tracks$interactions)$background.title <- 'black'
displayPars(tracks$gene)$background.title <- 'black'
displayPars(tracks$mut)$background.title <- 'black'
pdf(file="~/tmp/tnfaip3_ibd.pdf")
plotTracks(tracks,from=min(start(regions.gr)),to=138300000)
dev.off()
pdf(file="~/tmp/tnfaip3_ibd_full.pdf")
plotTracks(tracks)
dev.off()

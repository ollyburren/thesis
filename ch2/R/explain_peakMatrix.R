## figure to explain PCHi-C data

library(Gviz)
library(data.table)
library(biomaRt)
library(GenomicInteractions)

stuff<-fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full.txt")
bc <- stuff[,list(intcount=.N),by=baitID]
## for illustration use a bait with two separate interactions
##435074
sbait <- sample(bc[intcount==1,]$baitID,1)

pchic.DT <- stuff[baitID==sbait,]
e75.genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",  host="feb2014.archive.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
chr <- unique(pchic.DT$baitChr)
start <- min(c(pchic.DT$baitStart,pchic.DT$oeStart))
end <- max(c(pchic.DT$baitEnd,pchic.DT$oeEnd))

target<-BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=e75.genemart,
        chromosome=chr,
        start=start,
        end=end,
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.1, filters=list(with_ox_refseq_mrna=T))

cell.types <- names(pchic.DT)[11:28]
mdat <- melt(pchic.DT,id.vars=c("baitID","baitStart","baitEnd","oeID","oeChr","oeStart","oeEnd"),measure.vars=cell.types)[value>5,]

ideo<-fread("/rds/user/ob219/hpc-work/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))

tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=seqlevels(target),bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- target
ranges(tracks$gene)$features <- "protein_coding"
displayPars(tracks$gene)$fontsize.group <- 18
displayPars(tracks$gene)$protein_coding <- "black"
displayPars(tracks$gene)$fill<-"white"
displayPars(tracks$gene)$col<-"black"
#displayPars(tracks$gene)$stackHeight<-0.1
int.list <- split(mdat,as.character(mdat$variable))
#for(i in seq_along(int.list)){
cols <- c("black","firebrick1","dodgerblue")
col.cnt <- 1
for(i in sample(seq_along(int.list),3)){
  ct <- int.list[[i]]
  cell_type <- unique(as.character(ct$variable))
  bait.gr <- with(ct,GRanges(seqnames=Rle(paste("chr",oeChr,sep="")),IRanges(start=baitStart,end=baitEnd)))
  oe.gr <- with(ct,GRanges(seqnames=Rle(paste("chr",oeChr,sep="")),IRanges(start=oeStart,end=oeEnd)))
  interactions <- GenomicInteractions(bait.gr, oe.gr,counts=asinh(ct$value))
  tracks[[cell_type]] <- InteractionTrack(name=paste("Tissue",col.cnt), interactions)
  displayPars(tracks[[cell_type]])$anchor.height <- 0.3
  displayPars(tracks[[cell_type]])$col.anchors.line <- cols[col.cnt]
  displayPars(tracks[[cell_type]])$col.anchors.fill <- cols[col.cnt]
  displayPars(tracks[[cell_type]])$col.interactions <- cols[col.cnt]
  displayPars(tracks[[cell_type]])$background.title <- cols[col.cnt]
  col.cnt <- col.cnt + 1
}

pdf("~/tmp/pchic_explain.pdf")
plotTracks(tracks,from=start-7e4,to=end+7e4,sizes=c(1,2,3,3,3,3))
dev.off()

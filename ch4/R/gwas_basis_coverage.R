## basis stats GWAS

p.dt <- fread('/home/ob219/share/as_basis/GWAS/snp_manifest/gwas_june.tab')

library(cowplot)

ggplot(p.dt,aes(x=ref_a1.af)) + geom_histogram()
p.dt[,maf:=ifelse(ref_a1.af>0.5,1-ref_a1.af,ref_a1.af)]
ggplot(p.dt,aes(x=maf)) + geom_histogram()

library(GenomicRanges)
p.dt[,c('chr','position'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
p.gr <- with(p.dt,GRanges(seqnames=Rle(paste0("chr",chr)),ranges=IRanges(start=position,width=1L)))

library(karyoploteR)

pdf("~/tmp/cov_basis_1.pdf")
kp <- plotKaryotype(plot.type=1,chromosome=paste0("chr",c(seq(1,21,by=2),'X')))
ff <- kpPlotDensity(kp, data=p.gr,data.panel=1,col="firebrick1")
dev.off()

pdf("~/tmp/cov_basis_2.pdf")
kp <- plotKaryotype(plot.type=1,chromosome=paste0("chr",c(seq(2,22,by=2),'Y')))
ff <- kpPlotDensity(kp, data=p.gr,data.panel=1,col="firebrick1")
dev.off()

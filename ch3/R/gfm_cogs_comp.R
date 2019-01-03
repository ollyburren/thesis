## ichip summary stats plots

## how about we rerun pmi using the GFM regions so we get a proper comparison ?

library(data.table)
library(magrittr)

man <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/support/ichip_manifest.csv")
man[,c('samples','prop'):=list(cases+controls,cases/(cases+controls))]

rfiles <- list.files(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/support/dense.ic.regions.guessfm.shuffled',full.names=TRUE)

all.cmds <- lapply(rfiles,function(f){
lapply(split(man,man$trait),function(d){
  prop <- d$prop
  n_samples <- d$samples
  gwas_tbx <- file.path('/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/gwas',paste(d$filename,'gz',sep='.'))
  kg_dir <- '/home/ob219/rds/hpc-work/DATA/1kgenome/VCF/EUR/by.chr.phase3/ALL.'
  kg_suffix <- '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz'
  region_file <- f
  out_dir <- file.path('/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/',d$label,'/')
  if(!file.exists(out_dir))
    dir.create(out_dir)
  tabix_bin <- '/home/ob219/bin/htslib/tabix'
  RSCRIPT <- '/home/ob219/git/CHIGP/R/computePPi.R'
  cmd <- sprintf('Rscript %s --tabix_bin=%s --region_file=%s --out_dir=%s --gwas_tbx=%s  --gwas_type=CC --n_samples=%f --prop_cases=%f --kg_compress_dir=%s --kg_compress_suffix=%s --pi_i=0.0001 --do_pmi=1',
  RSCRIPT,tabix_bin,region_file,out_dir,gwas_tbx,n_samples,prop,kg_dir,kg_suffix)
}) %>% do.call('c',.)
}) %>% do.call('c',.)

write(all.cmds,'/home/ob219/tmp/rerun_ichip_pmi/cmds.txt')


all.cmds <- lapply(rfiles,function(f){
lapply(split(man,man$trait),function(d){
  prop <- d$prop
  n_samples <- d$samples
  gwas_tbx <- file.path('/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/gwas',paste(d$filename,'gz',sep='.'))
  kg_dir <- '/home/ob219/rds/hpc-work/DATA/1kgenome/VCF/EUR/by.chr.phase3/ALL.'
  kg_suffix <- '.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz'
  region_file <- f
  out_dir <- file.path('/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions/',d$label,'/')
  if(!file.exists(out_dir))
    dir.create(out_dir)
  tabix_bin <- '/home/ob219/bin/htslib/tabix'
  RSCRIPT <- '/home/ob219/git/CHIGP/R/computePPi.R'
  cmd <- sprintf('Rscript %s --tabix_bin=%s --region_file=%s --out_dir=%s --gwas_tbx=%s  --gwas_type=CC --n_samples=%f --prop_cases=%f --kg_compress_dir=%s --kg_compress_suffix=%s --pi_i=0.0001 --do_pmi=0',
  RSCRIPT,tabix_bin,region_file,out_dir,gwas_tbx,n_samples,prop,kg_dir,kg_suffix)
}) %>% do.call('c',.)
}) %>% do.call('c',.)

write(all.cmds,'/home/ob219/tmp/rerun_ichip_pmi/nopmi_cmds.txt')


stop("Run on q")

## compile into one file
out.dir <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed'
dirs <- list.dirs(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/raw',full.names=TRUE)[-1]
options("scipen"=100)
lapply(dirs,function(d){
  all.files <- list.files(path=d,pattern="*.pmi",full.names=TRUE)
  DT <- lapply(all.files,fread) %>% rbindlist
  DT <- DT[order(chr,start),]
  fname <- file.path(out.dir,paste(basename(d),'pmi',sep='.'))
  write.table(DT,file=fname,row.names=FALSE,sep="\t",quote=FALSE)

})
options("scipen"=0)

## compile into one file no pmi
out.dir <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions'
dirs <- list.dirs(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions/raw',full.names=TRUE)[-1]
options("scipen"=100)
lapply(dirs,function(d){
  all.files <- list.files(path=d,pattern="*.ppi",full.names=TRUE)
  DT <- lapply(all.files,fread) %>% rbindlist
  DT <- DT[order(chr,start),]
  fname <- file.path(out.dir,paste(basename(d),'ppi',sep='.'))
  write.table(DT,file=fname,row.names=FALSE,sep="\t",quote=FALSE)

})
options("scipen"=0)

## next get gene scores
pmi.files <- list.files(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed',pattern='*.pmi',full.names=TRUE)

cmds <- lapply(pmi.files,function(f){

pmi_file = f
out_dir = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/cogs'
int = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_interactions.RData'
frags = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_frags.by.ld.RData'
csnps = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_csnps.by.ld.RData'
sets = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/support/tnact_tree.yaml'


RSCRIPT <- '/home/ob219/git/CHIGP/R/computeGeneScoreH.R'
cmd <- sprintf('Rscript %s --pmi_file=%s --out_dir=%s --int=%s --frags=%s --csnps=%s target.gene.cSNPs.only=1 --sets=%s --include.interactions.only=1 --decompose=1 --ppi.thresh=0.01 --BF.thresh=3',
RSCRIPT,pmi_file,out_dir,int,frags,csnps,sets)
}) %>% do.call('c',.)
write(cmds,'/home/ob219/tmp/rerun_ichip_pmi/cogs.txt')

stop("Run")


pmi.files <- list.files(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions/processed',pattern='*.ppi',full.names=TRUE)

cmds <- lapply(pmi.files,function(f){

pmi_file = f
out_dir = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions/cogs'
int = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_interactions.RData'
frags = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_frags.by.ld.RData'
csnps = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/RDATA/javierre_tnact_gfm_regions_csnps.by.ld.RData'
sets = '/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/support/tnact_tree.yaml'


RSCRIPT <- '/home/ob219/git/CHIGP/R/computeGeneScoreH.R'
cmd <- sprintf('Rscript %s --pmi_file=%s --out_dir=%s --int=%s --frags=%s --csnps=%s target.gene.cSNPs.only=1 --sets=%s --include.interactions.only=1 --decompose=1 --ppi.thresh=0.01 --BF.thresh=3',
RSCRIPT,pmi_file,out_dir,int,frags,csnps,sets)
}) %>% do.call('c',.)
write(cmds,'/home/ob219/tmp/rerun_ichip_pmi/ppi_cogs.txt')

stop("Run")




library(data.table)
library(cowplot)
library(magrittr)

#DATA.DIR <- '/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi/'
DATA.DIR <- '/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed'

files <- list.files(path=DATA.DIR,pattern="*.pmi",full.names=TRUE)

ic.dat <- lapply(files,function(f){
  DT <- fread(f)
  DT[,disease:=gsub("\\_IC.pmi","",basename(f))]
}) %>% rbindlist

ic.dat[,list(count=.N),by=disease]

rbind(ic.dat[ppi<0.01 & is.na(imp.snp.pos),list(count=.N,cat='below'),by=disease],
ic.dat[ppi>0.01 & is.na(imp.snp.pos),list(count=.N,cat='above'),by=disease])


#ggplot(ic.dat,aes(x=disease)) + geom_bar()

#ABF.DIR <- '/Users/oliver/DATA/JAVIERRE_ICHIP/out/tnact_hierarchical_geneScore_0.01'
ABF.DIR <- '/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/cogs'

files <- list.files(path=ABF.DIR,pattern="*pmi_full.tab",full.names=TRUE)

ppi.thresh <- 0.01 # for table change to 0.5


abf <- lapply(files,fread) %>% rbindlist
abf[,disease:=sub("\\_IC","",disease)]
abf <- abf[disease %in% c('ATD','T1D','RA','CEL') & biotype=='protein_coding',]
abf <- abf[overall_gene_score>ppi.thresh,.(ensg=ensg,gscore=overall_gene_score,disease=disease,type='ABF')]

gfm <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/GUESSFM_JUNE_13/gene_prioritisation_0.01.csv")

#gfm <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/GUESSFM_FINAL/gene_prioritisation_0.01_csd3.csv")

gfm <- gfm[disease %in% c('ATD','T1D','RA','CEL') & overall_ppi>ppi.thresh,.(ensg=ensg,gscore=overall_ppi,disease=disease,type='GFM')]
setkeyv(gfm,c('ensg','disease'))
setkeyv(abf,c('ensg','disease'))
tot<-rbind(abf,gfm)
length(intersect(unique(tot[type=='ABF',]$ensg),unique(tot[type!='ABF',]$ensg)))
tot<-melt(tot,id.vars=c('ensg','disease','type'))
tot <- dcast(tot,ensg+disease~type+variable,fill=0)
tot[,label:=FALSE]
tot[ABF_gscore > ppi.thresh & GFM_gscore>ppi.thresh,label:=TRUE]

library(biomaRt)

ensembl_archive="feb2014.archive.ensembl.org" # Final V37 Ensembl
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
tot.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
filters="ensembl_gene_id",values=unique(tot$ensg), mart=ensembl) %>% data.table
setkey(tot.genes,ensembl_gene_id)
setkey(tot,ensg)
anno.tot <- tot[tot.genes]

pc.tot<-anno.tot[gene_biotype=='protein_coding',]




pp<-ggplot(pc.tot[ABF_gscore>0.5 | GFM_gscore>0.5],aes(x=ABF_gscore,y=GFM_gscore)) +
geom_point(size=1) + geom_abline(a=1,col='red',lty=2) + xlab("aBF COGS score") +
ylab("GUESSFM COGS score") + facet_wrap(~disease) + theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))


## ok use GUESSFM analysis to assign a gene to an LD block - will be rough and ready

load("/Users/oliver/DATA/JAVIERRE_ICHIP/out/GUESSFM_JUNE_13/MPPI-overall.RData")
agg <- melt(agg2,id.vars=c('ensg','ld.id'),measure.vars=c('ICOELIAC','RA','T1D','GRAVES'))
ld.group <- agg[agg[, .I[which.max(value)], by=c('ensg','variable')]$V1]
ld.group[variable=='ICOELIAC',variable:='CEL']
ld.group[variable=='GRAVES',variable:='ATD']
setkeyv(ld.group,c('ensg','variable'))
setkeyv(pc.tot,c('ensg','disease'))
pc.tot <- ld.group[pc.tot]
setnames(pc.tot,'variable','disease')

## assign gfm region based on proximity - gave up doing analytically

gfm.gr <- readRDS("/Users/oliver/DATA/JAVIERRE_ICHIP/support/GUESSFM.gr.RDS")
pgenes.gr <- with(pc.tot,GRanges(seqnames=Rle(chromosome_name),ranges=IRanges(start=start_position,end=end_position)))

pc.tot[,c('nearest.region','nearest.region.dist'):=list(gfm.gr$gfm_name[nearest(pgenes.gr,gfm.gr)],
  mcols(distanceToNearest(pgenes.gr,gfm.gr))$distance)
  ]

## what is the ranking of the different methods within disease ld block

pc.tot[,c('rGFM','rABF'):=list(rank(-GFM_gscore),rank(-ABF_gscore)),by=c('disease','nearest.region')]
# a summary of those where the top ranked gene agrees
pc.tot[,rankMatch:=FALSE]
pc.tot[rGFM<2 & rABF<2  & GFM_gscore!=0 & ABF_gscore!=0,rankMatch:=TRUE]

## add in where gfm results are missing as these are not comparable

pc.tot[,gfm.miss:=nearest.region.dist>2.5e6 & GFM_gscore==0]


## at 0.5 how many agree ?

table(pc.tot[ABF_gscore>0.5 | GFM_gscore>0.5,]$rankMatch)

athresh <- pc.tot[(ABF_gscore>0.5 | GFM_gscore>0.5) & nearest.region<3e6,]
table(athresh[,list(match=sum(rGFM<2 & rABF<2 & GFM_gscore!=0 & ABF_gscore!=0)),by=c('disease','nearest.region')]$match)
olsum<-athresh[,list(match=sum(rGFM<2 & rABF<2 & GFM_gscore!=0 & ABF_gscore!=0)),by=c('disease','nearest.region')]
olsum[,gene.region:=paste(disease,nearest.region,sep=':')]
pc.tot[,gene.region:=paste(disease,nearest.region,sep=':')]
pc.tot[,top.match.region:=FALSE]
pc.tot[gene.region %in% olsum[match!=0,]$gene.region,top.match.region:=TRUE]
pc.tot[,top.match:=FALSE]
pc.tot[,top.match:=rGFM<2 & rABF<2 & GFM_gscore!=0 & ABF_gscore!=0]
athresh <- pc.tot[ABF_gscore>0.5 | GFM_gscore>0.5,]



save_plot("~/tmp/cogs_comparison.pdf",pp)

athresh <- athresh[gfm.miss==FALSE,]

summ <- athresh[ensg %in% pc.tot$ensg & gfm.miss==FALSE,list(both=sum(ABF_gscore > 0.5 & GFM_gscore>0.5),
ABF=sum(ABF_gscore>0.5 & GFM_gscore<0.5),
GFM=sum(ABF_gscore<0.5 & GFM_gscore>0.5),tot=.N),by=disease]

library(xtable)
setcolorder(summ,c('disease','ABF','GFM','both','tot'))
xtable(summ[order(tot ,decreasing = TRUE),])

tot[,list(duplicated=sum(!duplicated(ensg))),by=disease]

## work out counts of genes per region per disease to create a stat

tp<-melt(athresh,id.vars=c('nearest.region','disease','ensg'),measure.vars=c('ABF_gscore','GFM_gscore'))[,list(count=sum(value>0.5)),by=c('nearest.region','disease','variable')]
tp[,label:=paste(disease,gsub("\\_gscore","",variable))]
tp[,group:=paste(nearest.region,disease)]

tp<-melt(tp,id.vars=c('disease','nearest.region','variable'),measure.vars='count') %>% dcast(.,disease+nearest.region~variable)
tp[,delta:=ABF_gscore-GFM_gscore]

pd <- position_dodge(0.4)

ggplot(tp,aes(x=disease,y=delta,color=disease)) + geom_jitter(height=0,width=0.1) + scale_y_continuous(breaks=seq(-8,8,by=1)) + guides(color=FALSE) +

## what is correlation between top ranking genes perhaps add to plot by color ?


## what is the correlation between two methods for both
both <- athresh[label==TRUE,]
# by disease

athresh[,list(cor=cor(GFM_gscore,ABF_gscore),p=cor.test(GFM_gscore,ABF_gscore)$p.value),by=disease]
both[,list(cor=cor(GFM_gscore,ABF_gscore),p=cor.test(GFM_gscore,ABF_gscore)$p.value),by=disease]


pp<-ggplot(pc.tot[(ABF_gscore>0.5 | GFM_gscore>0.5) & gfm.miss==FALSE, ],aes(x=ABF_gscore,y=GFM_gscore)) +
geom_point(size=1) + geom_abline(a=1,col='red',lty=2) + xlab("aBF COGS score") +
ylab("GUESSFM COGS score") + facet_wrap(~disease) + theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))

 #+ geom_hline(yintercept=0.5,col='steelblue',lty=2) +
#geom_vline(xintercept=0.5,col='steelblue',lty=2)

save_plot("~/tmp/cogs_comparison.pdf",pp)

## what about supers where aBF and GUESSFM agree on top gene and they overlap eRNA etc. !!!

prior <- fread("/Users/oliver/Downloads/13059_2017_1285_MOESM9_ESM.gz")
best<-athresh[label==TRUE & rankMatch==TRUE,]
prior[paste(Ensembl_GeneID,Disease) %in% paste(best$ensg,best$disease),topmatch:=TRUE]
takefor<-prior[!is.na(diff.expr) & Closest_Min_P_Value_Susceptibility_SNP_P_Value<5e-8 & Disease %in% c('T1D','CEL','RA','ATD') & analysis=='ICHIP_GF',.(Ensembl_GeneID,Gene_Name,Disease,diff.expr,diff.erna)] %>% unique
takefor<-takefor[order(Gene_Name),]
## add in gene position

tf.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
filters="ensembl_gene_id",values=unique(takefor$Ensembl_GeneID), mart=ensembl) %>% data.table
setkey(tf.genes,ensembl_gene_id)
setkey(takefor,Ensembl_GeneID)
takefor <- takefor[tf.genes]

tfgenes.gr <- with(takefor,GRanges(seqnames=Rle(chromosome_name),ranges=IRanges(start=start_position,end=end_position)))
## label closest region
takefor[,c('nearest.region','nearest.region.dist'):=list(gfm.gr$gfm_name[nearest(tfgenes.gr,gfm.gr)],
  mcols(distanceToNearest(tfgenes.gr,gfm.gr))$distance)
  ]


## get details

### NO IMPUTATION !!!!!!

# ABF.DIR <- '/Users/oliver/DATA/JAVIERRE_ICHIP/out/nopmi_with_gfm_regions/cogs'
#
# files <- list.files(path=ABF.DIR,pattern="*ppi_full.tab",full.names=TRUE)
#
# ppi.thresh <- 0.01 # for table change to 0.5
#
#
# abf <- lapply(files,fread) %>% rbindlist
# abf[,disease:=sub("\\_IC","",disease)]
# abf <- abf[disease %in% c('ATD','T1D','RA','CEL') & biotype=='protein_coding',]
# abf <- abf[overall_gene_score>ppi.thresh,.(ensg=ensg,gscore=overall_gene_score,disease=disease,type='ABF')]
#
# gfm <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/GUESSFM_JUNE_13/gene_prioritisation_0.01.csv")
#
# #gfm <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/GUESSFM_FINAL/gene_prioritisation_0.01_csd3.csv")
#
# gfm <- gfm[disease %in% c('ATD','T1D','RA','CEL') & overall_ppi>ppi.thresh,.(ensg=ensg,gscore=overall_ppi,disease=disease,type='GFM')]
# setkeyv(gfm,c('ensg','disease'))
# setkeyv(abf,c('ensg','disease'))
# tot<-rbind(abf,gfm)
# length(intersect(unique(tot[type=='ABF',]$ensg),unique(tot[type!='ABF',]$ensg)))
# tot<-melt(tot,id.vars=c('ensg','disease','type'))
# tot <- dcast(tot,ensg+disease~type+variable,fill=0)
# tot[,label:=FALSE]
# tot[ABF_gscore > ppi.thresh & GFM_gscore>ppi.thresh,label:=TRUE]
#
# ## some of these are not protein coding these need to be removed use biomart
#
#
#
#
# summ <- tot[,list(both=sum(ABF_gscore > ppi.thresh & GFM_gscore>ppi.thresh),
# ABF=sum(ABF_gscore>ppi.thresh & GFM_gscore<ppi.thresh),
# GFM=sum(ABF_gscore<ppi.thresh & GFM_gscore>ppi.thresh),tot=.N),by=disease]
#
# library(xtable)
# setcolorder(summ,c('disease','ABF','GFM','both','tot'))
# xtable(summ[order(tot ,decreasing = TRUE),])
#
# tot[,list(duplicated=sum(!duplicated(ensg))),by=disease]
#
#
# pp<-ggplot(tot[ABF_gscore>0.5 | GFM_gscore>0.5],aes(x=ABF_gscore,y=GFM_gscore)) +
# geom_point(size=1) + geom_abline(a=1,col='red',lty=2) + xlab("aBF COGS score") +
# ylab("GUESSFM COGS score") + facet_wrap(~disease) + theme(axis.text.x=element_text(angle = -90, hjust = 0, vjust=0.5))
#
#  #+ geom_hline(yintercept=0.5,col='steelblue',lty=2) +
# #geom_vline(xintercept=0.5,col='steelblue',lty=2)
#
# save_plot("~/tmp/cogs_comparison.pdf",pp)



# summ <- pc.tot[,list(both=sum(ABF_gscore > ppi.thresh & GFM_gscore>ppi.thresh),
# ABF=sum(ABF_gscore>ppi.thresh & GFM_gscore<ppi.thresh),
# GFM=sum(ABF_gscore<ppi.thresh & GFM_gscore>ppi.thresh),tot=.N),by=disease]
#
# library(xtable)
# setcolorder(summ,c('disease','ABF','GFM','both','tot'))
# xtable(summ[order(tot ,decreasing = TRUE),])




## lets investigate a case of where there is no overlap - use T1D

t1d <- pc.tot[disease=='T1D',]
t1d[,label:=FALSE]
t1d[ABF_gscore > 0.5 & GFM_gscore>0.5,label:=TRUE]
t1d.abf <- t1d[ABF_gscore>0.5 & label==FALSE,]



## chr19 region looks like a good one to focus on as it contains TYK2 a
## functionally validated gene
jav.DT <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab')
# reg.chr <- '19'
# reg.start <- 10.3e6
# reg.end <- 10628607
# setkey(t1d,ensg)
#
# rint <- t1d[chromosome_name==reg.chr & start_position>reg.start-100 & end_position<reg.end+100 & (ABF_gscore>0.01 | GFM_gscore>0.01) & gene_biotype=='protein_coding',  ]
#
# gfm.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=GFM_gscore))
# abf.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=ABF_gscore))
#
#
# ## load in marginals
#
# marg <- fread("~/Downloads/13059_2017_1285_MOESM7_ESM")
# t1d.marg <- marg[chr==reg.chr & between(position,reg.start,reg.end) & trait=='T1D',]
# marg.gr <- with(t1d.marg,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=position,width=1L),mppi=mppi))
#
# ## load in aBF
#
# abf.t1d <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi/T1D_IC.pmi")
# t1d.abf <- abf.t1d[chr==reg.chr & between(start,reg.start,reg.end),]
# abf.raw.gr <- with(t1d.abf,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start,width=1L),ppi=ppi))
#
# ## add in interaction tracks
#
#
# fjav.DT <- jav.DT[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,]
# fjav.DT <- fjav.DT[ensg %in% rint$ensg,]
# tbait.gr <- with(fjav.DT[ensg == 'ENSG00000105397',],GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd)))
#
#
# bait.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd),idx=1:nrow(fjav.DT)))
# oe.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=oeStart,end=oeEnd),idx=1:nrow(fjav.DT)))
# oe.idx <- subsetByOverlaps(oe.gr,tbait.gr)$idx
# interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=1)
# library(Gviz)
# library(GenomicInteractions)
# ideo<-fread("/Users/oliver/DATA/ucsc/hg19_cytoBandIdeo.txt")
# setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
# tracks<-list()
# tracks$itr <- IdeogramTrack(genome="hg19", chromosome=reg.chr,bands=ideo)
# tracks$axis<-GenomeAxisTrack()
# tracks$gene <- BiomartGeneRegionTrack(
#         genome="hg19", name="Genes", transcriptAnnotation="symbol",
#         mart=ensembl_archive,
#         gene = rint$ensg,
#         fill="black",
#         col="black",
#         collapseTranscripts="gene",shape = "arrow",
#         stackHeight=0.5, filters=list(ensembl_gene_id=rint$ensg,biotype='protein_coding'))
# tracks$interactions <- InteractionTrack(name='pcHi-C', interaction)
# displayPars(tracks$interactions)$anchor.height <- 0.1
# displayPars(tracks$interactions)$col.anchors.line <- "black"
# displayPars(tracks$interactions)$col.anchors.fill <- "darkgreen"
# displayPars(tracks$interactions)$col.interactions <- "darkgreen"
# displayPars(tracks$interactions)$col.outside <- "darkgreen"
#
# tracks$gfm.marg<- DataTrack(name="mPPi",marg.gr,col="steelblue")
# tracks$gfm<- DataTrack(name="GFM COGS",gfm.gr,type='histogram',col="steelblue",fill="steelblue")
# tracks$abf.ppi<- DataTrack(name="PPi",abf.raw.gr,col="firebrick")
# tracks$abf<- DataTrack(name="aBF COGS",abf.gr,type='histogram',col="firebrick",fill="firebrick")
#
# displayPars(tracks$axis)$background.title <- 'black'
# displayPars(tracks$gene)$background.title <- 'black'
# displayPars(tracks$gfm)$background.title <- 'steelblue'
# displayPars(tracks$gfm.marg)$background.title <- 'steelblue'
# displayPars(tracks$abf)$background.title <- 'firebrick'
# displayPars(tracks$abf.ppi)$background.title <- 'firebrick'
# displayPars(tracks$interactions)$background.title <- 'darkgreen'
#
# displayPars(tracks$gene)$fontsize <-20
# displayPars(tracks$gfm)$fontsize <-16
# displayPars(tracks$gfm.marg)$fontsize <-16
# displayPars(tracks$abf)$fontsize <-16
# displayPars(tracks$abf.ppi)$fontsize <-16
# displayPars(tracks$abf.ppi)$size <-50
# displayPars(tracks$gfm)$ylim <- c(0,1)
# displayPars(tracks$gfm.marg)$ylim <- c(0,1)
# displayPars(tracks$abf.ppi)$ylim <- c(0,1)
# displayPars(tracks$abf)$ylim <- c(0,1)
#
# pdf("~/tmp/abf_disc.pdf")
# plotTracks(tracks,sizes=c(0.5,0.5,1,0.5,rep(1,4)))
# dev.off()
#
#
# ## what about the other way around ?
#
# t1d.gfm <- t1d[GFM_gscore>0.5 & label==FALSE & gene_biotype=='protein_coding',]
# t1d.gfm <- t1d.gfm[order(chromosome_name,start_position),]
#
# #16p11.2 is missing from analysis in T1D why ? - Not in the LD file for aBF - check iChip region file against GUESSFM region files
# #chr16   28295305        29025978
#
# ## look at gfm and ichip and see which regions are missing and how much of the difference is explained
#
# abf.reg <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/support/dense.ic.regions.shuffled.bed") %>%
#   with(.,GRanges(seqnames=Rle(V1),ranges=IRanges(start=V2,end=V3)))
# gfm.reg <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/support/dense.ic.regions.guessfm.bed")%>%
#   with(.,GRanges(seqnames=Rle(V1),ranges=IRanges(start=V2,end=V3)))
#
# ol<-findOverlaps(abf.reg,gfm.reg) %>% as.matrix
#
# abf.reg[-ol[,1],]
# gfm.reg[-ol[,2],]



# reg.chr <- '16'
# reg.start <- 28467693
# reg.end <- 28550495
#
# abf.t1d <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed/T1D_IC.pmi")
# t1d.abf <- abf.t1d[chr==reg.chr & between(start,reg.start,reg.end),]




reg.chr <- '16'
reg.start <- 28467693
reg.end <- 28550495
setkey(t1d,ensg)
setkey(t1d.genes,ensembl_gene_id)

rint <- t1d[chromosome_name==reg.chr & start_position>reg.start-1e6 & end_position<reg.end+1e6 & (ABF_gscore>0.01 | GFM_gscore>0.01) & gene_biotype=='protein_coding',  ]

gfm.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=GFM_gscore))
abf.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=ABF_gscore))


## load in marginals

marg <- fread("~/Downloads/13059_2017_1285_MOESM7_ESM")
t1d.marg <- marg[chr==reg.chr & between(position,reg.start,reg.end) & trait=='T1D',]
marg.gr <- with(t1d.marg,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=position,width=1L),mppi=mppi))

## load in aBF

#abf.t1d <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi/T1D_IC.pmi")
abf.t1d <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed/T1D_IC.pmi")
t1d.abf <- abf.t1d[chr==reg.chr & between(start,reg.start,reg.end),]
abf.raw.gr <- with(t1d.abf,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start,width=1L),ppi=ppi))

## add in interaction tracks

fjav.DT <- jav.DT[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,]
fjav.DT <- fjav.DT[ensg %in% rint$ensg,]
tbait.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd)))
tbait.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd)))

bait.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd),idx=1:nrow(fjav.DT)))
oe.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=oeStart,end=oeEnd),idx=1:nrow(fjav.DT)))
oe.idx <- subsetByOverlaps(oe.gr,tbait.gr)$idx
interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=1)
library(Gviz)
library(GenomicInteractions)
ideo<-fread("/Users/oliver/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=reg.chr,bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=ensembl_archive,
        gene = rint$ensg,
        fill="black",
        col="black",
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.5, filters=list(ensembl_gene_id=rint$ensg,biotype='protein_coding'))
tracks$interactions <- InteractionTrack(name='PCHi-C', interaction)
displayPars(tracks$interactions)$anchor.height <- 0.1
displayPars(tracks$interactions)$col.anchors.line <- "black"
displayPars(tracks$interactions)$col.anchors.fill <- "darkgreen"
displayPars(tracks$interactions)$col.interactions <- "darkgreen"
displayPars(tracks$interactions)$col.outside <- "darkgreen"

tracks$gfm.marg<- DataTrack(name="mCVPP",marg.gr,col="steelblue")
tracks$gfm<- DataTrack(name="mCOGS",gfm.gr,type='histogram',col="steelblue",fill="steelblue")
tracks$abf.ppi<- DataTrack(name="sCVPP",abf.raw.gr,col="firebrick")
tracks$abf<- DataTrack(name="sCOGS",abf.gr,type='histogram',col="firebrick",fill="firebrick")

displayPars(tracks$axis)$background.title <- 'black'
displayPars(tracks$gene)$background.title <- 'black'
displayPars(tracks$gfm)$background.title <- 'steelblue'
displayPars(tracks$gfm.marg)$background.title <- 'steelblue'
displayPars(tracks$abf)$background.title <- 'firebrick'
displayPars(tracks$abf.ppi)$background.title <- 'firebrick'
displayPars(tracks$interactions)$background.title <- 'darkgreen'

displayPars(tracks$gene)$fontsize <-20
displayPars(tracks$gfm)$fontsize <-16
displayPars(tracks$gfm.marg)$fontsize <-16
displayPars(tracks$abf)$fontsize <-16
displayPars(tracks$abf.ppi)$fontsize <-16
displayPars(tracks$abf.ppi)$size <-50
displayPars(tracks$gfm)$ylim <- c(0,1)
displayPars(tracks$gfm.marg)$ylim <- c(0,1)
displayPars(tracks$abf.ppi)$ylim <- c(0,1)
displayPars(tracks$abf)$ylim <- c(0,1)

pdf("~/tmp/gfm_disc.pdf",useDingbats=FALSE)
plotTracks(tracks,sizes=c(0.5,0.5,1,0.5,rep(1,4)))
dev.off()

## what about 19 with non imputed data ?

t1d <- tot[disease=='T1D',]
t1d[,label:=FALSE]
t1d[ABF_gscore > 0.5 & GFM_gscore>0.5,label:=TRUE]
t1d.abf <- t1d[ABF_gscore>0.5 & label==FALSE,]

jav.DT <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full_denorm_bait2baits_e75.tab')
reg.chr <- '19'
reg.start <- 10.3e6
reg.end <- 10628607
setkey(t1d,ensg)

rint <- t1d[chromosome_name==reg.chr & start_position>reg.start-100 & end_position<reg.end+100 & (ABF_gscore>0.01 | GFM_gscore>0.01) & gene_biotype=='protein_coding',  ]

gfm.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=GFM_gscore))
abf.gr <- with(rint,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start_position,end=end_position),gfm=ABF_gscore))


## load in marginals

marg <- fread("~/Downloads/13059_2017_1285_MOESM7_ESM")
t1d.marg <- marg[chr==reg.chr & between(position,reg.start,reg.end) & trait=='T1D',]
marg.gr <- with(t1d.marg,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=position,width=1L),mppi=mppi))

## load in aBF

abf.t1d <- fread("/Users/oliver/DATA/JAVIERRE_ICHIP/out/pmi_with_gfm_regions/processed//T1D_IC.pmi")
t1d.abf <- abf.t1d[chr==reg.chr & between(start,reg.start,reg.end),]
abf.raw.gr <- with(t1d.abf,GRanges(seqnames=Rle(paste0('chr',reg.chr)),ranges=IRanges(start=start,width=1L),ppi=ppi))

## add in interaction tracks


fjav.DT <- jav.DT[Total_CD4_NonActivated>5 | Total_CD4_Activated>5,]
fjav.DT <- fjav.DT[ensg %in% rint$ensg,]
tbait.gr <- with(fjav.DT[ensg == 'ENSG00000105397',],GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd)))


bait.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=baitStart,end=baitEnd),idx=1:nrow(fjav.DT)))
oe.gr <- with(fjav.DT,GRanges(seqnames=Rle(paste0('chr',reg.chr)),IRanges(start=oeStart,end=oeEnd),idx=1:nrow(fjav.DT)))
oe.idx <- subsetByOverlaps(oe.gr,tbait.gr)$idx
interaction <- GenomicInteractions(bait.gr[oe.idx,], oe.gr[oe.idx,], counts=1)
library(Gviz)
library(GenomicInteractions)
ideo<-fread("/Users/oliver/DATA/ucsc/hg19_cytoBandIdeo.txt")
setnames(ideo,c('chrom', 'chromStart', 'chromEnd', 'name', 'gieStain'))
tracks<-list()
tracks$itr <- IdeogramTrack(genome="hg19", chromosome=reg.chr,bands=ideo)
tracks$axis<-GenomeAxisTrack()
tracks$gene <- BiomartGeneRegionTrack(
        genome="hg19", name="Genes", transcriptAnnotation="symbol",
        mart=ensembl_archive,
        gene = rint$ensg,
        fill="black",
        col="black",
        collapseTranscripts="gene",shape = "arrow",
        stackHeight=0.5, filters=list(ensembl_gene_id=rint$ensg,biotype='protein_coding'))
tracks$interactions <- InteractionTrack(name='PCHi-C', interaction)
displayPars(tracks$interactions)$anchor.height <- 0.1
displayPars(tracks$interactions)$col.anchors.line <- "black"
displayPars(tracks$interactions)$col.anchors.fill <- "darkgreen"
displayPars(tracks$interactions)$col.interactions <- "darkgreen"
displayPars(tracks$interactions)$col.outside <- "darkgreen"

tracks$gfm.marg<- DataTrack(name="mCVPP",marg.gr,col="steelblue")
tracks$gfm<- DataTrack(name="mCOGS",gfm.gr,type='histogram',col="steelblue",fill="steelblue")
tracks$abf.ppi<- DataTrack(name="sCVPP",abf.raw.gr,col="firebrick")
tracks$abf<- DataTrack(name="sCOGS",abf.gr,type='histogram',col="firebrick",fill="firebrick")

displayPars(tracks$axis)$background.title <- 'black'
displayPars(tracks$gene)$background.title <- 'black'
displayPars(tracks$gfm)$background.title <- 'steelblue'
displayPars(tracks$gfm.marg)$background.title <- 'steelblue'
displayPars(tracks$abf)$background.title <- 'firebrick'
displayPars(tracks$abf.ppi)$background.title <- 'firebrick'
displayPars(tracks$interactions)$background.title <- 'darkgreen'

displayPars(tracks$gene)$fontsize <-20
displayPars(tracks$gfm)$fontsize <-16
displayPars(tracks$gfm.marg)$fontsize <-16
displayPars(tracks$abf)$fontsize <-16
displayPars(tracks$abf.ppi)$fontsize <-16
displayPars(tracks$abf.ppi)$size <-50
displayPars(tracks$gfm)$ylim <- c(0,1)
displayPars(tracks$gfm.marg)$ylim <- c(0,1)
displayPars(tracks$abf.ppi)$ylim <- c(0,1)
displayPars(tracks$abf)$ylim <- c(0,1)

pdf("~/tmp/abf_disc.pdf",useDingbats=FALSE)
plotTracks(tracks,sizes=c(0.5,0.5,1,0.5,rep(1,4)))
dev.off()

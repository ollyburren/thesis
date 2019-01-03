RAW_DIR <- '/home/ob219/share/cogs_bb/summary_stats/raw/'
BB_PIDS <-

## create a master list of pid with associated reference MAFs and LD regions
if(FALSE){
  library(GenomicRanges)
  BB.DT <- fread("~/share/cogs_bb/bb_variant_list.txt")
  BB.DT[,c('chr','pos','a1','a2'):=tstrsplit(variant,':')]
  BB.DT <- BB.DT[chr %in% 1:22,]
  pids <- BB.DT[,pid:=paste(chr,pos,sep=':')]$pid
  ## assign to LD regions
  UK10K <- readRDS('/home/ob219/rds/hpc-work/DATA/UK10K/UK10K_0.005_MAF.RDS')
  UK10K[,pid:=paste(CHROM,POS,sep=':')]
  UK10K <- UK10K[pid %in% pids,]
  LD.DT <- fread("/home/ob219/rds/hpc-work/as_basis/support/all.1cM.bed",header=FALSE)[,c('chr','start'):=tstrsplit(V1,':')]
  LD.DT[,c('start','end'):=tstrsplit(start,'-')]
  LD.DT <- LD.DT[order(as.numeric(chr),as.numeric(start)),]
  ld.gr <- with(LD.DT,GRanges(seqnames=Rle(chr),ranges=IRanges(start=as.numeric(start),end=as.numeric(end))))
  ld.gr <- sort(ld.gr)
  uk.gr <- with(UK10K,GRanges(seqnames=Rle(CHROM),ranges=IRanges(POS,width=1L)))
  ol <- findOverlaps(uk.gr,ld.gr) %>% as.matrix
  UK10K[ol[,1],ld.block:=ol[,2]]
  ## remove MHC
  UK10K <- UK10K[!(CHROM=='6' & between(POS,25e6,35e6)),]
  UK10K <- UK10K[!duplicated(pid),.(pid,maf=MAF,ld.block)]
  saveRDS(UK10K,file="~/share/cogs_bb/snp_manifest.RDS")
  write(pids[!duplicated(pids)],"~/share/cogs_bb/bb_pid_list")
  ## assign
}else{
  UK10K <- readRDS("~/share/cogs_bb/snp_manifest.RDS")
}

## process SLE

SLE.DT <- file.path(RAW_DIR,'sle_benett_2016.tab') %>% fread
## keep analysis to non-autosomes only
SLE.DT <- SLE.DT[Chr<23,]
SLE.DT <- SLE.DT[,.(pid=paste(Chr,Distance,sep=':'),beta=Beta1Add,se.beta=Se1Add,p.value=PValueAdd)]
M <- merge(SLE.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/SLE.RDS")
rm(SLE.DT)

## process UC
UC.DT <- file.path(RAW_DIR,'uc_build37_45975_20161107.txt') %>% fread
UC.DT[,pid:=sub("\\_.*$","",MarkerName)]
UC.DT <- UC.DT[,.(pid,beta=Effect,se.beta=StdErr,p.value=P.value)]
M <- merge(UC.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/UC.RDS")
rm(UC.DT)

## process CD
CD.DT <- file.path(RAW_DIR,'cd_build37_40266_20161107.txt') %>% fread
CD.DT[,pid:=sub("\\_.*$","",MarkerName)]
CD.DT <- CD.DT[,.(pid,beta=Effect,se.beta=StdErr,p.value=P.value)]
M <- merge(CD.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/CD.RDS")
rm(CD.DT)

## process RA
RA.DT <- file.path(RAW_DIR,'RA_GWASmeta_European_v2.txt') %>% fread
setnames(RA.DT,names(RA.DT) %>% make.names)
## need to compute standard error of beta
RA.DT[,se.beta:=log(OR.A1.)/qnorm(P.val/2,lower.tail=FALSE)]
RA.DT[,beta:=log(OR.A1.)]
RA.DT <- RA.DT[,.(pid=paste(Chr,Position.hg19.,sep=':'),beta,se.beta,p.value=P.val)]
M <- merge(RA.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/RA.RDS")
rm(RA.DT)

## process PBC - cant use as not imputed to 1000 genomes

PBC.DT <- file.path(RAW_DIR,'PBC_GCMETA_fixedeffects') %>% fread
## on build 36
library(rtracklayer)
PBC.DT[,id:=1:.N]
pbc.36.gr <- with(PBC.DT,GRanges(seqnames=Rle(paste0('chr',chr)),ranges=IRanges(start=pos,width=1L),id=id))
c<-import.chain('/home/ob219/rds/hpc-work/DATA/LIFTOVER/hg18ToHg19.over.chain') ## e.g. hg19ToHg18.over.chain
pbc.37.gr<-unlist(liftOver(pbc.36.gr,c))
DT.37 <- data.table(id=pbc.37.gr$id,position.37=start(pbc.37.gr))
PBC.DT <- merge(PBC.DT,DT.37,by.x='id',by.y='id',all.x=TRUE)
## 173 snps don't match after coord conversion
PBC.DT <- PBC.DT[!is.na(position.37),]
PBC.DT <- PBC.DT[,.(pid=paste(chr,position.37,sep=':'),beta,se.beta=se,p.value=P_value)]
M <- merge(PBC.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/PBC.RDS")
rm(PBC.DT)

## process T1D

T1D.DT <- file.path(RAW_DIR,'t1d_cooper_2017.txt') %>% fread
T1D.DT <- T1D.DT[,.(pid=paste(chromosome,position,sep=':'),beta=beta.meta,se.beta=se.meta,p.value=p.meta)]
T1D.DT <- T1D.DT[!is.na(beta),]
M <- merge(T1D.DT,UK10K,by='pid')
M[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
M <- M[order(chr,pos),]
saveRDS(M,file="/home/ob219/share/cogs_bb/summary_stats/T1D.RDS")
rm(T1D.DT)

## could do clever BCF things but load all and get the overlap to compute a SNP snp_manifest
## we know that they are all in biobank as we already included that at the beginning

gf <- list.files(path="/home/ob219/share/cogs_bb/summary_stats",pattern="*.RDS",full.names=TRUE)

master_pid <- UK10K$pid
for(f in gf){
  message(f)
   pid <- readRDS(f)$pid
   message(sprintf("Length of disease pid %d",length(pid)))
   master_pid <- intersect(master_pid,pid)
   message(sprintf("Master PID length %d",length(master_pid)))
}

## using the above create a UK10K manifest file that can be used as a filter
UK10K.filter <- UK10K[pid %in% master_pid,]
## a useful addition is to match pids to UKBB variant identifiers
## grab a random bb summary stat
ukbb <- fread("zcat /home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/20001_1063.gwas.imputed_v3.both_sexes.tsv.bgz")
ukbb <- ukbb[,.(variant)][,pid:=sub("^([^:]+:[^:]+):.*","\\1",variant)]
ukbb <- ukbb[,.(bbid=variant,pid,order=1:.N)]
M <- merge(UK10K.filter,ukbb,by='pid')
UK10K.filter <- M[!duplicated(pid),]
UK10K.filter <- UK10K.filter[order(order),]
UK10K.filter <- UK10K.filter[,order:=NULL]
saveRDS(UK10K.filter,"/home/ob219/share/cogs_bb/snp_manifest.RDS")

man.DT <- readRDS("/home/ob219/share/cogs_bb/snp_manifest.RDS")
gf <- list.files(path="/home/ob219/share/cogs_bb/summary_stats",pattern="*.RDS",full.names=TRUE)
for(f in gf){
  message(f)
   DT <- readRDS(f)
   message(sprintf("Rows before filtering %d",nrow(DT)))
   DT<-DT[pid %in% man.DT$pid,]
   DT<-DT[!duplicated(pid),]
   message(sprintf("Rows after filtering %d",nrow(DT)))
   saveRDS(DT,file=f)
}

## create svCPP for each trait that have the same format as the

logsum <- function(x) {
  my.max <- max(x) ##take out the maximum value in log form)
  my.res <- my.max + log(sum(exp(x - my.max )))
  return(my.res)
}

wakefield_pp <- function(p,f, N, s,pi_i=1e-4,sd.prior=0.2) {
    if(length(p) != length(f))
      stop("p and f must be vectors of the same size")
    # compute V
    V <- 1 / (2 * N * f * (1 - f) * s * (1 - s))
    # convert p vals to z
    z <- qnorm(0.5 * p, lower.tail = FALSE)
    ## Shrinkage factor: ratio of the prior variance to the total variance
    r <- sd.prior^2 / (sd.prior^2 + V)
    ## Approximate BF
    lABF = 0.5 * (log(1-r) + (r * z^2))
    ## tABF - to add one we create another element at the end of zero for which pi_i is 1
    tABF <- c(lABF,0)
    vpi_i<-c(rep(pi_i,length(lABF)),1)
    sBF <- logsum(tABF + log(vpi_i))
    exp(lABF+log(pi_i)-sBF)
}

trait.meta <- list(
  CD.RDS = list(n1=12194,n0=28072),
  RA.RDS = list(n1=14361,n0=43923),
  SLE.RDS = list(n1=4036,n0=6959),
  T1D.RDS = list(n1=5913,n0=8829),
  UC.RDS = list(n1=12366,n0=33609)
)

gf <- list.files(path="/home/ob219/share/cogs_bb/summary_stats",pattern="*.RDS",full.names=TRUE)
for(f in gf){
  DT <- readRDS(f)
  ## check for pvalues = 0 or 1
  DT[p.value==1,p.value:=0.99]
  DT[p.value==0,p.value:=1e-200]
  ## next compute
  meta <- trait.meta[[basename(f)]]
  N <- meta$n1 + meta$n0
  s <-  meta$n1/N
  DT[,cvPP:=wakefield_pp(p.value,maf,N,s),by=ld.block]
  message(summary(DT$cvPP))
  saveRDS(DT[,.(chr,pos,cvPP,ld.block)],file=file.path('/home/ob219/share/cogs_bb/sCVPP/',basename(f)))
  message(sprintf("Wrote %s",file.path('/home/ob219/share/cogs_bb/sCVPP/',basename(f))))
}

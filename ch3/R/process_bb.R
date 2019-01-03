#

OUT_DIR <- '/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/cogs_bb_posteriors/'

if(FALSE){
  ## code to create a bb trait manifest with stuff linked to codes
  bb_dir <- '/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018'
  bb_files <- list.files(path=bb_dir,pattern="*.tsv.bgz",full.names=TRUE)
  bb.DT <- data.table(share_file=bb_files,share_name=basename(bb_files))
  bb_phenofile<-'/home/ob219/rds/hpc-work/as_basis/bb/bb_gwas_link_list.20180731.csv'
  pheno <- fread(bb_phenofile)
  setnames(pheno,names(pheno) %>% make.names)
  pheno <- pheno[Phenotype.Code!='N/A',.(File,Phenotype.Code,Phenotype.Description)]
  pheno[,phe:=make.names(Phenotype.Description) %>% gsub("Non.cancer.illness.code..self.reported..|Treatment.medication.code..|Cancer.code..self.reported..","",.)]
  M <- merge(bb.DT,pheno,by.y='File',by.x='share_name')
  ## load in phenotype file
  P <- fread('/home/ob219/rds/hpc-work/as_basis/bb/summary_stats_20180731/phenotypes.both_sexes.tsv')
  P<-P[,.(phenotype,non_missing=n_non_missing,cases=n_cases,controls=n_controls)]
  M <- merge(M,P,by.y='phenotype',by.x='Phenotype.Code')
  saveRDS(M,"/home/ob219/share/cogs_bb/bb_trait_manifest.RDS")
}

if(FALSE){
  bb.trait.DT <- readRDS("/home/ob219/share/cogs_bb/bb_trait_manifest.RDS")
  SCRIPT <- '/home/ob219/git/thesis/process_bb.R'
  if(file.exists(out.file))
  cmd <- sapply(bb.trait.DT$share_file,function(f){
    p <- bb.trait.DT[share_file==args$fname]
    out.file <- file.path(OUT_DIR,paste(p$phe,'RDS',sep='.'))
    if(!file.exists(out.file))
      sprintf("Rscript %s -f %s",SCRIPT,f)
  })
  write(cmd,file="~/tmp/qstuff/bb_cvpp.txt")
}

library(optparse)

TEST<-FALSE
option_list = list(
        make_option(c("-f", "--fname"), type="character", default=NULL,
              help="index of phenotype to process ", metavar="character")
            )
if(!TEST){
  opt_parser = OptionParser(option_list=option_list);
  args = parse_args(opt_parser)
  if (is.null(args$fname)){
	   print_help(opt_parser)
	    stop("Supply a vcf file to process", call.=FALSE)
    }
}else{
  args <- list(fname='/home/ob219/share/Data/GWAS-summary/uk_biobank_neale_summary_stats_2018/20003_2038460150.gwas.imputed_v3.both_sexes.tsv.bgz')
}

bb.trait.DT <- readRDS("/home/ob219/share/cogs_bb/bb_trait_manifest.RDS")
p <- bb.trait.DT[share_file==args$fname]
out.file <- file.path(OUT_DIR,paste(p$phe,'RDS',sep='.'))
message(out.file)
if(file.exists(out.file))
  q(save="no")
man.DT <- readRDS("/home/ob219/share/cogs_bb/snp_manifest.RDS")
DT <- sprintf("zcat %s",args$fname) %>% fread
DT <- DT[variant %in% man.DT$bbid,]

computeOR <- function(n,n1,Sx,Sxy) {
    ## estimated allele freqs in cases and controls
    fe1 <- Sxy/(2*n1)
    fe0 <- (Sx - Sxy)/(2*(n-n1))
    ## estimated odds ratio
    fe1 * (1-fe0) / ( (1-fe1) * fe0 )
}

SElor<-function(n,n1,Sx,Sxy){
    n0<-n-n1
    #fe1 is the af in cases
    c <- Sxy/(2*n1)
    #fe0 is af in the controls
    a <- (Sx - Sxy)/(2*(n-n1))
    b<-1-a
    d<-1-c
    ## normalise
    a<-(a*n0)/n
    b<-(b*n0)/n
    c<-(c*n1)/n
    d<-(d*n1)/n
    ## estimated odds ratio bc/ad
    sqrt(1/2) * sqrt(1/n) * sqrt(1/a + 1/b + 1/c + 1/d)
}

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


p <- bb.trait.DT[share_file==args$fname]
DT[,or:=computeOR(p$non_missing,p$cases,AC,ytx)]
DT[,c('theta','se.theta'):=list(log(or),SElor(p$non_missing,p$cases,AC,ytx))]
DT[,c('theta.pval','theta.Z','n0','n1'):=list(2*(pnorm(abs(theta/se.theta),lower.tail = FALSE)),theta/se.theta,p$non_missing-p$cases,p$cases)]
## note that the OR method for deriving p.values has numerical issues use quant p.vals instead
out <- DT[,.(pid=sub("^([^:]+:[^:]+):.*","\\1",variant),beta=theta,se.beta=se.theta,p.value=pval)]
out <- merge(out,man.DT[,.(pid,maf,ld.block)],by='pid')
out[,c('chr','pos'):=tstrsplit(pid,':') %>% lapply(.,as.numeric)]
out <- out[order(chr,pos)]
## finally add causal variant posterior probabilities.
N <- p$non_missing
s <- p$cases/N
out[,cvPP:=wakefield_pp(p.value,maf,N,s),by=ld.block]
saveRDS(out[,.(chr,pos,cvPP,ld.block)],file=out.file)
message(sprintf("Wrote %s",out.file))

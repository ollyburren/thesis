### replot stahl and okada figure and examine distribution of imputed and unimputed values
## for this analysis we take variants from chromosome 1 and select those variants that are
## typed in stahl et al. and select these variants from okada - and redo PMI.

## we then have the pmi imputed values and the actual pmi values that we can compare.

##step 1 obtain a list of variants to filter okada on.
stahl.DT <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/out/pmi/RA.pmi")
## unimputed (by PMI) are imp.snp.pos is na
s.DT<-stahl.DT[is.na(imp.snp.pos) & chr==1,]
s.DT[,pid:=paste(chr,start,sep=':')]
## read in okada dataset and filter based on the above

okada.DT <- fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/gwas/OKADA_RA.bed")
setnames(okada.DT,c('chr','start','end','id','p.value'))
o.DT <- okada.DT[chr==1,]
o.DT[,pid:=paste(chr,start,sep=':')]
o.DT <- o.DT[pid %in% s.DT$pid,]

write.table(o.DT,file="/home/ob219/rds/hpc-work/thesis/ra_pmi_analysis/ra_okada_stah_chr1.bed",col.names=FALSE,row.names=FALSE,quote=FALSE,sep="\t")

## tried to run on interactive node and q
#Rscript --vanilla PMI.R --region_file=~/rds/hpc-work/thesis/ra_pmi_analysis/0.1_chr1.txt --out_dir=/home/ob219/rds/hpc-work/thesis/ra_pmi_analysis/out --gwas_tbx=/home/ob219/rds/hpc-work/thesis/ra_pmi_analysis/ra_okada_stah_chr1.bed.gz --gwas_type=CC --n_samples=58284 --prop_cases=0.25 --kg_compress_dir=/home/ob219/rds/hpc-work/DATA/1kgenome/VCF/EUR/by.chr.phase3/ALL. --kg_compress_suffix=.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.EUR.vcf.gz --tabix_bin=/home/ob219/bin/htslib/tabix --pi_1=0.0001

## load in these results and compare to the results from okada where they did 'proper' imputation

pmi.DT <- fread("/home/ob219/rds/hpc-work/thesis/ra_pmi_analysis/out0.1_chr1.txt.pmi")
pmi.DT[,pid:=paste(chr,start,sep=':')]

pmi.DT<-pmi.DT[,.(pid,pmi.p=-log10(pval),pmi.imputed=!is.na(imp.r2))]
omi.DT <- o.DT[,.(pid,okada.mlp=-log10(p.value))]
omi.DT[,source:=pid %in% s.DT$pid]


M<-merge(pmi.DT,omi.DT,by='pid')

## is the distribution of imputed SNPs in okada different from those that are non imputed ?
library(cowplot)
ggplot(M,aes(x=qnorm((10^-okada.mlp)/2,lower.tail=FALSE),color=source)) + geom_density() + coord_cartesian(xlim=c(0,5))


ppa <- ggplot(M[pmi.imputed==TRUE,],aes(x=okada.mlp,y=pmi.p)) +
geom_hex(bins=60) +
ylab(expression("PMI"~-log[10](P))) +
xlab(expression("Okada et al."~-log[10](P))) +
geom_abline(a=1,alpha=0.5,lty=2,color='red')
save_plot("~/tmp/ppi_comp.pdf",ppa)


ppb <- ggplot(M[pmi.imputed==TRUE,],aes(x=qnorm((10^-okada.mlp)/2,lower.tail=FALSE),y=qnorm((10^-pmi.p)/2,lower.tail=FALSE))) +
#geom_point(alpha=0.5) +

ylab(expression("PMI"~-log[10](P))) +
xlab(expression("Okada et al."~-log[10](P))) +
geom_abline(a=1,alpha=0.5,lty=2,color='red')

plot_grid(ppa,ppb)

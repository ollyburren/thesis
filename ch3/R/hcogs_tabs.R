library(data.table)
library(magrittr)

files <- list.files(path="/Users/oliver/DATA/THESIS/ch3/hgenescore_0.01",pattern="*pmi_prioritised.tab",full.names=TRUE)
files <- files[-grep("COOPER_T1D.pmi_prioritised.tab",files)]

all.dat <- lapply(files,fread) %>% rbindlist

## filter so only protein coding genes

GTF<-fread('/Users/oliver/DATA/JAVIERRE_GWAS/support/Homo_sapiens.GRCh37.75.genes.gtf')
GTF[,ensg:=gsub('.*"(.*)"',"\\1",tstrsplit(GTF$V9,';',fixed=TRUE,fill=TRUE)[[1]])]
# only protein coding
all.dat<-all.dat[ensg %in% GTF[V2=='protein_coding',]$ensg,]

## create a summary table for each disease and label

## perhaps only look in AI disease ?

aid <- c('CD_IMB','CEL','MS','PBC','RA_OKADA_IMB','SLE','T1D','UC')

sum<-all.dat[disease %in% aid & isLeaf==TRUE ,list(cat.count=.N),by=c('node','disease')]
sum <- melt(sum,id.vars=c('node','disease'),measure.vars='cat.count')
mat <- dcast(sum,node~disease,fill=0)

mat.mat <- as.matrix(mat[,-1])
rownames(mat.mat) <- mat$node
#mat.mat <- log(mat.mat + 1)
library(xtable)
xtable(mat.mat)

## what are coding genes prioritised ?

all.dat[disease %in% aid & node=='coding',]

all.dat[disease %in% aid,]

## need a script that takes a list of genes and a tissue and returns the fragment level scores -
in this way we can see what is driving the category assignments and also see how many are due to multiple genes being
prioritised by a single fragment.



## check out the virtual promoter assignments.
hind3 <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/support/HindIII_baits_e75.bed')
setnames(hind3,c('chr','start','end','id','deet'))
hind3[,c('gname','ensg','biotype','strand'):=tstrsplit(deet,':')]
hind3.gene<-hind3[biotype=='protein_coding',]

hind3 <- fread('/Users/oliver/DATA/JAVIERRE_GWAS/support/Digest_Human_HindIII.bed')

## what is the average size of promoter coding regions ?

## get a list of promoter prioritised genes

p.genes <- unique(all.dat[disease %in% aid & node == 'promoter',]$ensg)

prom.id<-hind3.gene[ensg %in% p.genes,]$id
prom.id <- split(c(prom.id-1,prom.id,prom.id+1),prom.id)
hind3[,width:=V3-V2+1]
lapply(prom.id,function(x){
  hind3[V4 %in% x,]$width %>% sum
}) %>% do.call('c',.) %>% mean

## what is going on with UC ?

p.genes.uc <- unique(all.dat[disease =='UC' & node == 'promoter',]$ensg)


library(biomaRt)

ensembl_archive="feb2014.archive.ensembl.org" # Final V37 Ensembl
ensembl <- useMart(host=ensembl_archive, biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
e75.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),filters="ensembl_gene_id",values=p.genes.uc, mart=ensembl) %>% data.table
e75.genes[order(chromosome_name,start_position),]

## they don't appear to all cluster together.

## what about neut findings ?

p.genes.neut <- all.dat[disease %in% aid & node == 'Neutrophils',]$ensg
neut.genes<-getBM(c("ensembl_gene_id","gene_biotype","external_gene_id","hgnc_symbol", "chromosome_name","start_position", "end_position","strand"),
filters="ensembl_gene_id",values=p.genes.neut, mart=ensembl) %>% data.table
neut.genes[order(chromosome_name,start_position),]

setkey(neut.genes,ensembl_gene_id)
dg <- all.dat[disease %in% aid & node == 'Neutrophils',.(disease,ensg)]
setkey(dg,ensg)

dg <- dg[neut.genes]
lapply(split(dg,dg$disease),function(x){
  x[order(chromosome_name,start_position),]
})

## load in the full sets so that we can lookup the scores for the selected node.



library(cowplot)

#ggplot(all.dat,aes(x=node,y=appi)) + geom_boxplot() + theme(axis.text.x=element_text(angle = -90, hjust = 0))

files <- list.files(path="/Users/oliver/DATA/THESIS/ch3/hgenescore_0",pattern="*pmi_full.tab",full.names=TRUE)
files <- files[-grep("COOPER_T1D.pmi_full.tab",files)]

full.dat <- lapply(files,fread) %>% rbindlist


m.dt <- melt(full.dat[biotype=="protein_coding" & disease %in% aid,],id.vars=c("disease","ensg","name","biotype","strand","baitChr"))
## need to change category so we can join

levels(m.dt$variable)<-gsub("_prey_only_gene_score","",levels(m.dt$variable))
levels(m.dt$variable)<-gsub("set\\.","",levels(m.dt$variable))
levels(m.dt$variable)<-gsub("_gene_score","",levels(m.dt$variable))

m.dt[,variable:=as.character(variable)]

setkeyv(m.dt,c('disease','variable','ensg'))
setkeyv(all.dat,c('disease','node','ensg'))

with.score<-m.dt[all.dat][!is.na(value),]
bymean.ai<-with.score[,list(mscore=mean(value)),by=variable]
#with.score[,fv:=factor(variable,levels=bymean[order(mscore),]$variable)]
bymean.ai[,variable:=factor(variable,levels=bymean.ai[order(mscore),]$variable)]
bymean.ai[,ai:=TRUE]
#ppai <- ggplot(with.score[isLeaf==TRUE,],aes(x=fv,y=value)) + geom_boxplot() +
#theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
#xlab("COGS Category") + ylab("COGS Category Score") +
#geom_hline(yintercept=0.5,col='firebrick')



#ppai <- ggplot(bymean[variable %in% with.score[isLeaf==TRUE & !variable %in% c('coding','promoter'),]$variable,],aes(x=variable,y=mscore)) + geom_bar(stat="identity") + theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5))


## what about the others ?

## randomly select 8 traits

all.diseases <- unique(full.dat$disease)

rand.traits <- sample(all.diseases[!all.diseases %in% aid],8)

m.dt <- melt(full.dat[biotype=="protein_coding" & disease %in% rand.traits,],id.vars=c("disease","ensg","name","biotype","strand","baitChr"))
## need to change category so we can join

levels(m.dt$variable)<-gsub("_prey_only_gene_score","",levels(m.dt$variable))
levels(m.dt$variable)<-gsub("set\\.","",levels(m.dt$variable))
levels(m.dt$variable)<-gsub("_gene_score","",levels(m.dt$variable))

m.dt[,variable:=as.character(variable)]

setkeyv(m.dt,c('disease','variable','ensg'))
setkeyv(all.dat,c('disease','node','ensg'))

with.score<-m.dt[all.dat][!is.na(value),]
bymean.nonai<-with.score[,list(mscore=median(value)),by=variable]
with.score[,fv:=factor(variable,levels=bymean[order(mscore),]$variable)]
bymean.nonai[,variable:=factor(variable,levels=bymean.ai[order(mscore),]$variable)]
bymean.nonai[,ai:=FALSE]


all <- rbind(bymean.ai,bymean.nonai)
vlabs <- c(Endothelial_precursors='EndP',
Macrophages_M1=expression(M*phi*1),
Macrophages_M0=expression(M*phi*0),
Macrophages_M2=expression(M*phi*2),
Erythroblasts='Ery',
Neutrophils='Neut',
Total_CD8='tCD8',
Foetal_thymus='FetT',
Total_CD4_NonActivated='naCD4',
Megakaryocytes='MK',
Naive_CD4='nCD8',
Total_CD4_MF='tCD4',
Naive_B='nB',
Total_B='tB',
Naive_CD4='nCD4',
Monocytes='Mon',
Naive_CD8 = 'nCD8',
Total_CD4_Activated='aCD4')


toplot <- all[variable %in% with.score[isLeaf==TRUE & !variable %in% c('coding','promoter'),]$variable,]
#toplot[,variable:=factor(variable,levels=toplot[1:17,]$variable,labels=vlabs)]


ppi <- ggplot(toplot,aes(x=variable,y=mscore,fill=ai)) +
geom_bar(stat="identity",position="dodge") +
scale_fill_manual(guide=FALSE,values=c('steelblue','firebrick1')) +
theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
ylab("Mean COGS Score") + xlab("COGS Tissue Category") +
scale_x_discrete(labels = vlabs)
save_plot("~/tmp/hscore.pdf",ppi)

## bootstrap

all.diseases <- unique(full.dat$disease)
all.diseases <- all.diseases[all.diseases!='RA']
bs.dt <- melt(full.dat[biotype=="protein_coding",],id.vars=c("disease","ensg","name","biotype","strand","baitChr"))
levels(bs.dt$variable)<-gsub("_prey_only_gene_score","",levels(bs.dt$variable))
levels(bs.dt$variable)<-gsub("set\\.","",levels(bs.dt$variable))
levels(bs.dt$variable)<-gsub("_gene_score","",levels(bs.dt$variable))
setkeyv(bs.dt,c('disease','variable','ensg'))
setkeyv(all.dat,c('disease','node','ensg'))
bs.dt<-bs.dt[all.dat][!is.na(value),]

bs.res <- lapply(1:1000,function(i){
  rand.traits <- sample(all.diseases,8)
  bs.dt[disease %in% rand.traits,list(mscore=mean(value)),by=variable]
}) %>% rbindlist

bsmean.nonai <- bs.res[,list(mscore=mean(mscore),varscore=var(mscore)),by=variable]
bsmean.nonai[,variable:=factor(variable,levels=bymean.ai[order(mscore),]$variable)]
bsmean.nonai[,ai:=FALSE]

all <- rbind(bymean.ai,bsmean.nonai,fill=TRUE)

toplot <- all[variable %in% with.score[isLeaf==TRUE & !variable %in% c('coding','promoter'),]$variable,]

toplot <- all[variable %in% with.score[isLeaf==TRUE,]$variable,]


ppi <- ggplot(toplot,aes(x=variable,y=mscore,fill=ai)) +
geom_bar(stat="identity",position="dodge") +
scale_fill_manual(guide=FALSE,values=c('steelblue','firebrick1')) +
theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
ylab("Mean COGS Score") + xlab("COGS Tissue Category") +
scale_x_discrete(labels = vlabs)


## perhaps rerun with the correct tree structure (i.e. the one that we generate from PCA analysis)

OUT_DIR <- '/home/ob219/rds/hpc-work/thesis/hgenescore'
INTERACTIONS_FILE <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/RDATA/javierre_interactions.RData'
FRAG_FILE <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/RDATA/javierre_frags.by.ld.RData'
CSNP_FILE <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/RDATA/javierre_csnps.by.ld.RData'
SET_FILE <- '/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/support/thesis_tree.yaml'

cmd <- 'Rscript /home/ob219/git/CHIGP//R/computeGeneScoreH.R --pmi_file=%s --out_dir=%s --int=%s --frags=%s --csnps=%s --target.gene.cSNPs.only=1 --sets=%s --include.interactions.only=1 --decompose=1 --ppi.thresh=0.5 --BF.thresh=3'
pmi.files <- list.files(path='/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/out/pmi',pattern="*.pmi",full.names=TRUE)
cmds <- sapply(pmi.files,function(p){
  sprintf(cmd,p,OUT_DIR,INTERACTIONS_FILE,FRAG_FILE,CSNP_FILE,SET_FILE)
})

write(cmds,file="~/tmp/thesis_h_genescore.txt")


## where is PMI DIR ?

Rscript /home/ob219/git/CHIGP//R/computeGeneScoreH.R --pmi_file=/home/ob219/git//cd4chic/gwas_paper/DATA//out/pmi/RA_OKADA_IMB.pmi --out_dir=/home/ob219/git//cd4chic/gwas_paper/DATA//out/hierarchical_geneScore_tnact_0.01/ --int=/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_tnact_interactions.RData --frags=/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_tnact_frags.by.ld.RData --csnps=/home/ob219/scratch/DATA/JAVIERRE_GWAS/RDATA/javierre_tnact_csnps.by.ld.RData --target.gene.cSNPs.only=1 --sets=/home/ob219/scratch/DATA/JAVIERRE_ICHIP/support/tnact_tree.yaml --include.interactions.only=1 --decompose=1 --ppi.thresh=0.01 --BF.thresh=3

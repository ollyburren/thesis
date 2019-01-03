library(data.table)
library(magrittr)
library(cowplot)

files <- list.files(path="/Users/oliver/DATA/THESIS/ch3/hgenescore_0",pattern="*pmi_full.tab",full.names=TRUE)
files <- files[-grep("COOPER_T1D.pmi_full.tab",files)]

full.dat <- lapply(files,fread) %>% rbindlist

files <- list.files(path="/Users/oliver/DATA/THESIS/ch3/hgenescore_0",pattern="*pmi_prioritised.tab",full.names=TRUE)
files <- files[-grep("COOPER_T1D.pmi_prioritised.tab",files)]

all.dat <- lapply(files,fread) %>% rbindlist

all.diseases <- unique(full.dat$disease)
all.diseases <- all.diseases[all.diseases!='RA']
bs.dt <- melt(full.dat[biotype=="protein_coding",],id.vars=c("disease","ensg","name","biotype","strand","baitChr"))
levels(bs.dt$variable)<-gsub("_prey_only_gene_score","",levels(bs.dt$variable))
levels(bs.dt$variable)<-gsub("set\\.","",levels(bs.dt$variable))
levels(bs.dt$variable)<-gsub("_gene_score","",levels(bs.dt$variable))
setkeyv(bs.dt,c('disease','variable','ensg'))
setkeyv(all.dat,c('disease','node','ensg'))
bs.dt<-bs.dt[all.dat][!is.na(value),]

## do AI

aid <- c('CD_IMB','CEL','MS','PBC','RA_OKADA_IMB','SLE','T1D','UC')

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

bs.res <- lapply(1:10000,function(i){
  rand.traits <- sample(all.diseases,8)
  bs.dt[disease %in% rand.traits,list(mscore=mean(value)),by=variable]
}) %>% rbindlist

bsmean.nonai <- bs.res[,list(mscore=mean(mscore),varscore=var(mscore)),by=variable]
bsmean.nonai[,variable:=factor(variable,levels=bymean.ai[order(mscore),]$variable)]
bsmean.nonai[,ai:=FALSE]

all <- rbind(bymean.ai,bsmean.nonai,fill=TRUE)

toplot <- all[variable %in% with.score[isLeaf==TRUE & !variable %in% c('coding','promoter'),]$variable,]


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
Total_CD4_Activated='aCD4',
coding='Coding',
promoter='VProm')


ggplot(toplot,aes(x=variable,y=mscore,fill=ai)) +
geom_bar(stat="identity",position="dodge") +
scale_fill_manual(guide=FALSE,values=c('steelblue','firebrick1')) +
theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
ylab("Mean COGS Score") + xlab("COGS Tissue Category") +
scale_x_discrete(labels = vlabs)


setkey(bymean.ai,variable)
setkey(bsmean.nonai,variable)
res <- bymean.ai[bsmean.nonai]
res[,z:=(mscore-i.mscore)/sqrt(varscore)]
res[,p:=pnorm(abs(z),lower.tail=FALSE)]
#resplot <- res[variable %in% with.score[isLeaf==TRUE & !variable %in% c('coding','promoter'),]$variable,.(variable,z,p)]
resplot <- res[variable %in% with.score[isLeaf==TRUE,]$variable,.(variable,z,p)]

resplot[,variable:=factor(variable,levels=resplot[order(z),]$variable)]
resplot[,pv:=-log10(p) * sign(z)]

sig.thresh <- qnorm(0.05/nrow(resplot),lower.tail=FALSE)

ppi<-ggplot(resplot,aes(x=variable,y=z)) + geom_bar(stat="identity",col="black",fill="white") +
theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) +
ylab("Z Score") + xlab("hCOGS Tissue Category") +
scale_x_discrete(labels = vlabs) + geom_hline(yintercept=sig.thresh,col="firebrick1",lty=2)
save_plot("~/tmp/bootstrap_hscore.pdf",ppi)

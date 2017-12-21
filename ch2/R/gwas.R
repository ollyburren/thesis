library(data.table)
library(cowplot)

## CODE TO GENERATE SUMMARY STATS ON GWAS DATASETS USED IN JAVIERRE ET AL


DATA.DIR<-'/Users/oliver/DATA/JAVIERRE_GWAS/gwas/'
all.gwas<-list.files(path=DATA.DIR,pattern='*.bed.gz$',full.names=TRUE)


## to process each data set we want to count input (SNPs - i.e. not imputed that are above MAF threshold and in 1KG)

processFile<-function(f){
    message(f)
    cmd<-sprintf("gunzip -c %s",f)
    nrow(fread(cmd))
}

det<-lapply(all.gwas,processFile)

names(det)<-sub("(.*).bed.gz","\\1",basename(all.gwas))


## next get the number of SNPs that come out of PMI.

PMI.DIR<-'/Users/oliver/DATA/JAVIERRE_GWAS/out/pmi/'

all.pmi<-list.files(path=PMI.DIR,pattern='*.pmi$',full.names=TRUE)


processPMI<-function(f){
    message(f)
    DT<-fread(f)
    trait<-sub("(.*)\\.pmi$","\\1",basename(f))
    n.non<-sum(is.na(DT$imp.snp.pos))
    n.imp<-sum(!is.na(DT$imp.snp.pos))
    total<-n.non+n.imp
    data.table(trait=trait,n=n.non,n.imputed=n.imp,total=total)
}

det2<-lapply(all.pmi,processPMI)
det2<-rbindlist(det2)

## get current table that we want to update

tab<-fread("/Users/oliver/DATA/JAVIERRE_GWAS/support/gwas_manifest.csv")
## add reference
setkey(tab,label)
setkey(det2,trait)

ftab<-tab[det2]
ftab<-ftab[ftab$label!='RA',]
ftab<-ftab[order(ftab$category,ftab$n)]
ftab[,perc.imp:=round((n.imputed/total)*100)]
ftab[ftab$category=='T2D',]$category<-'Other'
ftab[ftab$label=='CD_IMB',label:='CRO']
ftab[ftab$label=='RA_OKADA_IMB',label:='RA']

library(reshape2)

bl<-melt(ftab,id.vars=c('label','category'),measure.vars=c('n','n.imputed'))
bl$label<-factor(gsub("\\_"," ",bl$label),levels=gsub("\\_"," ",ftab$label))
bl$variable<-factor(bl$variable,levels=c('n.imputed','n'))

lab.col<-list(Autoimmune='dodgerblue',Blood='firebrick',Metabolic='seagreen',Other='black')
a<-unlist(lab.col[bl$category])

library(ggplot2)
library(cowplot)
#pdf("~/git/thesis/ch2/pdf/pmi_gwas.pdf")

pp <- ggplot(bl[label!='COOPER T1D',],aes(x=label,y=value,fill=variable)) + geom_bar(color='black',stat='identity') +
    scale_fill_manual(name = 'PMI imputed',label=c('Yes','No'),values=c('grey','white')) +
    theme(axis.text.x=element_text(angle = -45, hjust = 0.0, vjust=1)) +
        theme(axis.text.x=element_text(colour = a)) + ylab("# SNPs") + xlab("Trait") + geom_hline(yintercept =median(ftab$total),color='red')
save_plot("~/git/thesis/ch2/pdf/pmi_gwas.pdf",pp,base_height=8)
#dev.off()

stop()
de.tab<-ftab[,.(category,trait,label,cases,controls,pmid)]
library('rentrez')
pdat<-lapply(unique(de.tab$pmid),function(p){
    tmp<-entrez_summary(db="pubmed", id=p)
    sprintf("%s et al. (%s)",tmp$authors[1,1],strsplit(tmp$pubdate,' ')[[1]][1])
})
names(pdat)<-unique(de.tab$pmid)
de.tab$author<-pdat[as.character(de.tab$pmid)]
library(xtable)
xtable(de.tab)

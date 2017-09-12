## joy plot
DATA.DIR<-'/Users/oliver/DATA/JAVIERRE_GWAS/out/geneScore'
fs<-list.files(path=DATA.DIR,pattern="*.tab",full.names = TRUE)
all.dat<-rbindlist(lapply(fs,fread))
all.dat.pc<-all.dat[all.dat$biotype=='protein_coding',]

anno<-fread('/Users/oliver/DATA/JAVIERRE_GWAS/support/gwas_manifest.csv')
anno<-subset(anno,label != 'RA')
anno[anno$label=='T2D',]$category='Other'
anno[anno$label=='CD_IMB',]$label<-'CRO'
anno[anno$label=='RA_OKADA_IMB',]$label<-'RA'


mean(sapply(split(all.dat.pc,all.dat.pc$disease),nrow))

library(ggplot2)
library(cowplot)

hs<-all.dat.pc[all.dat.pc$all_gene_score>0.5 & disease !='RA',]
hs[hs$disease=='RA_OKADA_IMB',]$disease<-'RA'
hs[hs$disease=='CD_IMB',]$disease<-'CRO'
setkey(anno,label)
setkey(hs,disease)
tp<-anno[hs]
counts<-tp[,list(gc=.N),by=c('label','category','cases','controls')]
counts<-counts[order(counts$category,counts$gc,decreasing = TRUE),]
counts$sample.size<-ifelse(!is.na(counts$controls),counts$cases + counts$controls,counts$cases)
tp$label<-factor(tp$label,levels=counts$label)
pdf("~/git/thesis/ch2/pdf/cogs_overall_summary.pdf")
ggplot(tp,aes(x=label,fill=category)) + geom_bar()  + scale_fill_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position=c(0.80, 0.80), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white")) + ylab('# Protein Coding Genes Prioritised') + xlab('Trait')
dev.off()

ggplot(counts,aes(x=log10(sample.size),y=gc)) + geom_point(aes(colour=category)) + xlab("Log10 (Sample Size)") + ylab ('# Prioritised genes') + geom_smooth(method='lm',formula = y ~ x,se=FALSE) + scale_color_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) + theme(axis.text.x=element_text(angle = -90, hjust = 0)) + theme(legend.position=c(0.1, 0.80), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))

quantile(all.dat.pc$all_gene_score,probs=seq(0,1,0.01))

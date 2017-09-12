#pdf(plot.file,width=9,height=7)
library(ggplot2)
library(data.table)
library(cowplot)
m.results<-readRDS("/Users/oliver/DATA/THESIS/ch2/bs_ch2.RDS")
m.results[m.results$lm.gwas=='T2D',]$category='Other'
m.results[m.results$lm.gwas=='CD',]$lm.gwas<-'CRO'
m.results[!m.results$category %in% c('Autoimmune','Blood'),]$lm.gwas<-''

re<-data.table(ymin=c(-3,0),ymax=c(5.5,5.5),xmin=c(0,-2.5),xmax=c(7,7),category=factor(c(1,2)))
#re<-re[1,]
pdf("~/git/thesis/ch2/pdf/bs_scatter_lm_vs_mmme.pdf")
p<-ggplot(re,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category)) + geom_rect(alpha=0.2) + scale_fill_manual(name = NULL,values =c('blue','green'),guide = FALSE)

p + geom_point(data=m.results,aes(x=lm.z,y=mmme.z ,pch=category,colour=category,label=lm.gwas),size=3,inherit.aes=FALSE) +
    theme_bw() +  xlab("< Myeloid vs Lymphoid Z Score >") +
    #ylab("<  (MON, MAC & NEU) vs (MEG & ERY) Z Score >") +
    ylab("<  Innate Immune vs Megakaryocyte/Erythroblast Z Score >") +
    geom_text(data=m.results,aes(x=lm.z,y=mmme.z ,pch=category,colour=category,label=lm.gwas),size=5,hjust=-0.2, vjust=0, angle=0, show_guide=FALSE,inherit.aes = FALSE) +
    scale_colour_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) + scale_shape_discrete(name=NULL,guide=FALSE) + coord_flip() +# ylim(c(-3,5.5)) +
    geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) + guides(size=FALSE,text=FALSE) +
    theme(legend.position=c(0.80, 0.80), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))
dev.off()

## compute wilcoxon for for Z scores of ai vs the rest on lymphoid

computeWilcoxon<-function(cat){
    message(cat)
    m.results$wilcox.cat<-as.character(m.results$category)
    m.results[m.results$wilcox.cat!=cat,]$wilcox.cat<-'Foo'
    m.results$wilcox.cat<-factor(m.results$wilcox.cat)
    
    data.table(lm=p.adjust(wilcox.test(data=m.results,lm.z ~ wilcox.cat,alternative='t')$p.value,n=8),
    mmme=p.adjust(wilcox.test(data=m.results,mmme.z ~ wilcox.cat,alternative='t')$p.value,n=8),cat=cat)
}

rbindlist(lapply(c('Autoimmune','Blood','Metabolic','Other'), computeWilcoxon))



#pdf(plot.file,width=9,height=7)
library(ggplot2)
library(data.table)
library(cowplot)
m.results<-readRDS("/Users/oliver/DATA/THESIS/ch2/bs_ch2.RDS")
m.results[m.results$lm.gwas=='T2D',]$category='Other'
m.results[!m.results$category %in% c('Autoimmune','Blood'),]$lm.gwas<-''

re<-data.table(ymin=c(-3,0),ymax=c(5.5,5.5),xmin=c(0,-2.5),xmax=c(7,7),category=factor(c(1,2)))
#re<-re[1,]

p<-ggplot(re,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category)) + geom_rect(alpha=0.2) + scale_fill_manual(name = NULL,values =c('blue','green'),guide = FALSE)

p + geom_point(data=m.results,aes(x=lm.z,y=mmme.z ,pch=category,colour=category,label=lm.gwas),size=3,inherit.aes=FALSE) +
    theme_bw() +  xlab("< Myeloid vs Lymphoid Z Score >") +
    ylab("<  (MON, MAC & NEU) vs (MEG & ERY) Z Score >") +
    geom_text(data=m.results,aes(x=lm.z,y=mmme.z ,pch=category,colour=category,label=lm.gwas),size=5,hjust=-0.2, vjust=0, angle=0, show_guide=FALSE,inherit.aes = FALSE) +
    scale_colour_manual(name = NULL,values = c(Autoimmune="red",Blood="black",Metabolic="darkgreen",Other="blue",T2D="yellow")) + scale_shape_discrete(name=NULL) + coord_flip() +# ylim(c(-3,5.5)) +
    geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) + guides(size=FALSE,text=FALSE) +
    theme(legend.position=c(0.85, 0.85), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))
    
expression(symbol('\256')
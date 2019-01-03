#pdf(plot.file,width=9,height=7)
library(ggplot2)
library(data.table)
library(cowplot)
library(ggrepel)
m.results<-readRDS("/Users/oliver/DATA/THESIS/ch2/bs_ch2.RDS")
m.results<-readRDS("/home/ob219/tmp/bs_ch2.RDS") %>% data.table
m.results[m.results$lm.gwas=='T2D',]$category='Other'
m.results[m.results$lm.gwas=='CD',]$lm.gwas<-'CRO'
<<<<<<< HEAD

re<-data.table(ymin=c(-3,0),ymax=c(5.5,5.5),xmin=c(0,-2.5),xmax=c(7,7),category=factor(c(1,2)))
p<-ggplot(re,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category)) + geom_rect(alpha=0.1) + scale_fill_manual(name = NULL,values =c('blue','green'),guide = FALSE)
p <- p + geom_point(data=m.results,aes(x=lm.z,y=mmme.z ,colour=category),size=3,inherit.aes=FALSE) +
    xlab("< Myeloid vs Lymphoid Z Score >") +
    ylab("<  Innate Immune vs Megakaryocyte/Erythroblast Z Score >") +
    scale_colour_manual(name = "Trait Category",values = c(Autoimmune="dodgerblue",Blood="firebrick",Metabolic="seagreen",Other="black")) +
    scale_shape_discrete(name=NULL,guide=FALSE) + coord_flip() +
    geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) +
    geom_text_repel(data=m.results,aes(x=lm.z,y=mmme.z,colour=category,label=lm.gwas),size=5,inherit.aes = FALSE,show.legend = FALSE) +
    background_grid(major = "xy", minor = "none") + scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
    theme(legend.position=c(0.70, 0.80), legend.text=element_text(size=15),legend.background=element_rect(fill = "white"),legend.key=element_rect(colour='black'))
p
save_plot("/Users/oliver/git/thesis/ch2/pdf/bs_scatter_lm_vs_mmme.pdf",p,base_height=7)
=======
#m.results[!m.results$category %in% c('Autoimmune','Blood'),]$lm.gwas<-''
m.results[!category %in% c('Autoimmune','Blood'),lm.gwas:="",]
re<-data.table(ymin=c(-5,0),ymax=c(3,3),xmin=c(0,-4),xmax=c(7,7),category=factor(c(1,2)))
#re<-re[1,]
#pdf("~/git/thesis/ch2/pdf/bs_scatter_lm_vs_mmme.pdf")
p<-ggplot(re,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category)) + geom_rect(alpha=0.2) + scale_fill_manual(name = NULL,values =c('blue','green'),guide = FALSE)
ppa <- p + geom_point(data=m.results,aes(x=lm.z,y=mmme.z * -1 ,colour=category),size=3,inherit.aes=FALSE) +
    theme_bw() +  xlab("< Myeloid vs Lymphoid Z Score >") +
    #ylab("<  (MON, MAC & NEU) vs (MEG & ERY) Z Score >") +
    ylab("< Erythroblast/Megakaryocyte vs Innate Immune Z Score >") +
    geom_text_repel(data=m.results,aes(x=lm.z,y=mmme.z * -1,colour=category,label=lm.gwas),size=5, show_guide=FALSE,inherit.aes = FALSE) +
    scale_colour_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) + scale_shape_discrete(name=NULL,guide=FALSE) + coord_flip() +# ylim(c(-3,5.5)) +
    geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) + guides(size=FALSE,text=FALSE) +
    theme(legend.position=c(0.20, 0.70), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))
#dev.off()
>>>>>>> b8a9719b0c5a3345320d124ab6c09fac2edcd1aa

#, legend.background = element_rect(colour = "white"), legend.key = element_rect(colour = "white")) +
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
<<<<<<< HEAD
=======


## regenerate the plot from burren et al.

#gwas x.z y.z cat max.z
library(data.table)
library(magrittr)
files <- list.files(path='/home/ob219/rds/hpc-work/burren_data/out/blockshifter_gap_5/',pattern="*.txt",full.names=TRUE)
dat <- lapply(files,fread) %>% rbindlist
dat[test=='Tnact',test:='X']
dat[test=='Total_CD4_NonActivated',test:='Y']
dat <- dat[test!="Control",]
dat <- melt(dat,id.vars=c('gwas','test'),measure.vars='z')
dat <- dcast(dat,gwas~test+variable)
dat <- dat[gwas!="RA",]
## cheat a bit here and use previous to get categories etc.

dat <- cbind(dat,m.results[,.(lm.gwas,category)])[,.(lm.gwas,X_z,Y_z,category)]
dat[lm.gwas=='T2D',]$category='Other'
dat[lm.gwas=='CD',]$lm.gwas<-'CRO'

re<-data.table(ymin=c(-2,0),ymax=c(5.5,5.5),xmin=c(0,-2.6),xmax=c(5.5,5.5),category=factor(c(1,2)))
#re<-re[1,]
#pdf("~/git/thesis/ch2/pdf/bs_scatter_lm_vs_mmme.pdf")
p<-ggplot(re,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=category)) + geom_rect(alpha=0.2) + scale_fill_manual(name = NULL,values =c('blue','green'),guide = FALSE)
ppb <- p + geom_point(data=dat,aes(x=X_z,y=Y_z * -1 ,colour=category),size=3,inherit.aes=FALSE) +
    theme_bw() +  xlab("< CD4+ Non Act vs CD4+ Act Z Score >") +
    #ylab("<  (MON, MAC & NEU) vs (MEG & ERY) Z Score >") +
    ylab("< Megakaryocyte/Erythroblast vs CD4+ Act/Nact Z Score >") +
    geom_text_repel(data=dat,aes(x=X_z,y=Y_z * -1,colour=category,label=lm.gwas),size=5, show_guide=FALSE,inherit.aes = FALSE) +
    scale_colour_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) + scale_shape_discrete(name=NULL,guide=FALSE) + coord_flip() +# ylim(c(-3,5.5)) +
    geom_vline(xintercept=0,alpha=0.3) + geom_hline(yintercept=0,alpha=0.3) + guides(size=FALSE,text=FALSE) +
    theme(legend.position=c(0.2, 0.7), legend.background = element_rect(colour = "grey"),legend.text = element_text(size = 15), legend.key = element_rect(colour = "white"))
#dev.off()

pdf("~/git/thesis/ch2/pdf/bs_scatter.pdf")
p <- plot_grid(ppa,ppb,labels = "AUTO", align = 'h')
save_plot("~/git/thesis/ch2/pdf/bs_scatter.pdf",p,base_width=10,base_height=5)
dev.off()

computeWilcoxon<-function(cat){
    message(cat)
    dat$wilcox.cat<-as.character(dat$category)
    dat[wilcox.cat!=cat,]$wilcox.cat<-'Foo'
    dat$wilcox.cat<-factor(dat$wilcox.cat)
    data.table(X=p.adjust(wilcox.test(data=dat,X_z ~ wilcox.cat,alternative='t')$p.value,n=8),
    Y=p.adjust(wilcox.test(data=dat,Y_z ~ wilcox.cat,alternative='t')$p.value,n=8),cat=cat)
}

rbindlist(lapply(c('Autoimmune','Blood','Metabolic','Other'), computeWilcoxon))


## addition of cooper !

all.files<-list.files(path='/home/ob219/scratch/burren_data/out/blockshifter_gap_5_panc',pattern="*.txt",full.names=TRUE)

library(magrittr)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggrepel)

res <- lapply(all.files,fread) %>% rbindlist
res <- res[,.(gwas,test,control,z)]
res[,comb:=paste(test,control,sep=':')]
res <- res[comb!='Control:PancIs',]
res<-melt(res[comb=='Total_CD4_NonActivated:Total_CD4_Activated',comb:='nact:act'],id.vars=c('gwas','comb'),measure.vars='z')
pl <- dcast(res,gwas~comb)
setnames(pl,make.names(names(pl)))

#ggplot(pl,aes(x=Control.PancIs,y=nact.act,label=gwas)) + geom_point() + geom_text_repel() + coord_flip()
pl[,gwas:=factor(gwas,levels=gwas[order(ME.PancIs)])]
#color this in based on trait type crib off the stuff above
tmpdata<-readRDS("~/tmp/bs_ch2.RDS")  %>% data.table
tmpdata<-tmpdata[,.(lm.gwas,category)]
setkey(tmpdata,lm.gwas)
setkey(pl,gwas)

pl<-tmpdata[pl]

pl[lm.gwas %in% c('GLUCOSE','GLUCOSE_BMI','INS_BMI','T2D'),category:='Metabolic']
pl[lm.gwas %in% c('COOPER_T1D','RA_OKADA_IMB','CD_IMB'),category:='Autoimmune']
pl[is.na(category),category:='Other']
pl[,lm.gwas:=factor(lm.gwas,levels=lm.gwas[order(ME.PancIs)])]

ggplot(pl[!lm.gwas %in% c('T1D','RA')],aes(x=lm.gwas,y=ME.PancIs,fill=category)) + geom_bar(stat="identity") +
theme(axis.text.x=element_text(angle = -90, hjust = 0,vjust=0.5)) + ylab("<PancIs Vs Ery,Mega Z>") +
scale_fill_manual(name = "Trait Category",values = c(Autoimmune="steelblue",Blood="firebrick",Metabolic="goldenrod2",Other="seagreen3")) + xlab("Trait")
>>>>>>> b8a9719b0c5a3345320d124ab6c09fac2edcd1aa

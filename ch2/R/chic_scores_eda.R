library(data.table)
library(magrittr)
stuff<-fread("/home/ob219/rds/hpc-work/DATA/JAVIERRE_GWAS/chic/merged_samples_12Apr2015_full.txt")

cell.types <- names(stuff)[12:28]

#M <- stuff[sample.int(.N,1e5),cell.types,with=FALSE] %>% as.matrix %>% dist %>% hclust


## compare distribution or arcsinh vs raw chicago scores

mDT <- melt(stuff,id.vars='oeID',measure.vars=cell.types)
setnames(mDT,c('oeID','cell_type','raw.chicago.score'))
mDT[,asinh.chicago.score:=asinh(raw.chicago.score)]

library(ggplot2)
library(ggridges)
library(cowplot)

#ppa <- ggplot(mDT[cell_type=="Total_B",],aes(x=pmin(raw.chicago.score,50))) + geom_histogram(bins=4000,col="black",fill="black") +
#coord_cartesian(xlim=c(0,50)) + xlab("CHiCAGO score") + ylab("PIR count") + geom_vline(xintercept=5,col='red',lty=2)

#ppa <- ggplot(mDT,aes(x=raw.chicago.score)) + geom_histogram(bins=4000,col="black",fill="black") +
#coord_cartesian(xlim=c(0,50)) + xlab("CHiCAGO score") + ylab("PIR count") + geom_vline(xintercept=5,col='red',lty=2)

ppa <- ggplot(mDT,aes(x=pmin(raw.chicago.score,50))) + geom_histogram(bins=100,col="black",fill="black") +
xlab("CHiCAGO score") + ylab("PIR count") + geom_vline(xintercept=5,col='firebrick1',lty=2)

ppb <- ggplot(mDT,aes(y=gsub("_"," ",cell_type),x=asinh.chicago.score)) + geom_density_ridges() +
geom_vline(xintercept=asinh(5),col="firebrick1",lty=2) + xlab("arsinh(CHiCAGO Score)") + ylab("Tissue Type")

ppc <- plot_grid(ppa,ppb,labels=c('A','B'),ncol=1)
save_plot("~/tmp/chic.score.pdf",ppc)
save_plot("~/tmp/chic.score.pdf",ppc,base_height=7)

prcomp <- stuff[,cell.types,with=FALSE] %>% as.matrix %>% apply(.,1,asinh) %>% prcomp


install.packages("ggdendro")

par(mfrow=c(1,1))
mat <- prcomp$x
rownames(mat) <- c('Mon','M0','M1','M2','Neu','MK','EndP','Ery','FetT','nCD4','tCD4','aCD4','naCD4','nCD8','tCD8','nB','tB')

library(ggdendro)
comp <- dist(mat) %>% hclust(.,method="complete") %>% as.dendrogram %>% dendro_data(., type = "rectangle")

comp$labels$y <- min(comp$segments$y) * 0.85
comp$segments$yend[comp$segments$yend==0] <- min(comp$segments$y) * 0.9

p <- ggplot(segment(comp)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_text(data = comp$labels, aes(x = x, y = y, label = label), size = 4,angle=90,hjust="right") + ylab("Height") +
  theme(axis.title.x=element_blank(),
       axis.text.x=element_blank(),
       axis.ticks.x=element_blank(),
     axis.line.x = element_blank()) + coord_cartesian(ylim=c(max(comp$segments$y),min(comp$segments$y) * 0.65))



save_plot("~/tmp/pchic_arcsinh_dend.pdf",p)
  #coord_flip() +
  #scale_y_reverse(expand = c(0.2, 0))

# dist(mat) %>% hclust(.,method="complete") %>% plot
# dist(prcomp$x) %>% hclust(.,method="average") %>% plot
# dist(prcomp$x) %>% hclust(.,method="ward.D") %>% plot
#
# ## what happens if we set anything below 5 to 0 ?
#
# test<- stuff[,cell.types,with=FALSE] %>% as.matrix
# #test[test<5]<-0
# #test <- pmin(test,5)
#
# test[test>2]<-0
# 
#
#
# prcomp2 <- apply(test,1,asinh) %>% prcomp
#
# dist(prcomp2$x) %>% hclust %>% plot

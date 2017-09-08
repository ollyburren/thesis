library(ggplot2)
library(data.table)

## this is code to generate Type I error simulations for blockshifter
## when scratch is back up get this.
(load("/scratch/ob219/share/blockshifter_reports.RData"))
all.sim.reports$sim<-factor(all.sim.reports$sim,levels=c('null','alternative','alternative2'))
all.sim.reports[all.sim.reports$platform=='_FULL',]$platform<-'NONE'
all.sim.reports$platform<-sub("\\_","",all.sim.reports$platform)
all.sim.reports$platform<-factor(all.sim.reports$platform,levels=c('NONE','HAPMAP','IL550'))

ft.dat<-data.table(platform=c('NONE','HAPMAP','IL550'),val=c(c(0.31425,0.38525,0.44475)))
ft.dat$sim <- "null"

df <- subset(all.sim.reports,tc=='combined')
df <- merge(df,ft.dat,by=c("platform","sim"),all=TRUE)
#df[sim!="null",val:=NA]
table(df$sim)
df$sim <- factor(as.character(df$sim),levels=c("alternative2","alternative","null")) # so facets in right order
df$platform <- factor(as.character(df$platform),levels=c("NONE","HAPMAP","IL550")) # so facets in right order
levels(df$platform) <- c("1,000 Genomes","HapMap","Illumina 550k")


cbbPalette <- c("#E69F00", "#56B4E9", "#009E73","#F0E442",
                "#0072B2", "#D55E00", "#CC79A7")


pdf("../pdf/blockshifter_type1_sims.pdf")
ggplot(df,aes(x=gap.size,y=obsd,ymin=ci.low,ymax=ci.hi,color=sim,shape=sim,group=sim)) +
geom_hline(aes(yintercept=val,color=sim),linetype=5) +
geom_linerange(size=0.5) +
geom_point(size=2) +
geom_line() +
theme_bw()  +
scale_colour_manual(name = "Scenario",values=cbbPalette[c(3,5,6)]) +
scale_linetype_discrete(guide=FALSE) +
scale_shape_manual(name="Simulation",values=c(15,17,19),guide=FALSE) +
xlab("SuperBlock Gap Size") + ylab("Type 1 error rate / power") +
facet_grid(.~platform,scales="free") +
geom_hline(yintercept = 0.05,linetype=5,colour="grey") +
theme(legend.position=c(0.9,0.9))
dev.off()

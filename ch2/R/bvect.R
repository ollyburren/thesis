library(cowplot)
library(data.table)
library(magrittr)

labs <- c('Coding','VProm','nB','tB','FetT','aCD4','naCD4','tCD4','nCD8','nCD4','tCD8','Mon','Neu','M2','M1','M0','EndP','MK','Ery')

all <- rep(TRUE,length(labs))
dt <- data.table(ct=factor(labs,levels=labs),score=all,instance=factor('Overall',levels=c('naCD4','Lymphoid','Overall')))

one <- logical(length(labs))
one[7] <- TRUE
dt.one <- data.table(ct=factor(labs,levels=labs),score=one,instance=factor('naCD4',levels=c('naCD4','Lymphoid','Overall')))

lymph <- one
lymph[3:11] <- TRUE

dt.lymph <- data.table(ct=factor(labs,levels=labs),score=lymph,instance=factor('Lymphoid',levels=c('naCD4','Lymphoid','Overall')))
all.DT <- rbindlist(list(dt,dt.one,dt.lymph))
all.DT[,label:=ifelse(score==TRUE,"1","0")]
col <- rep('black',length(labs))
col[1:2] <- 'grey'

ppb<-ggplot(all.DT,aes(y=instance,x=ct,fill=score,label=label)) + geom_tile(col='black') +
scale_fill_manual(name="Present",values=c('TRUE'='grey','FALSE'='white'),guide=FALSE) +
theme(axis.title.x = element_text(size=14),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1,colour=col),line=element_blank()) +
xlab("pcHi-C Tissue") + geom_text(col="black") +
ylab("'b' vector")

save_plot("~/tmp/bvect.pdf",ppb,base_height=2,base_width=7)

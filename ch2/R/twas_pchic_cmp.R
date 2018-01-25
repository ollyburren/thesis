library(data.table)
library(GenomicRanges)

## code to generate plaform stats for
foo<-list(coloc = "ENSG00000103811", twas = c("ENSG00000005020",
"ENSG00000196502", "ENSG00000110917", "ENSG00000137261", "ENSG00000089127",
"ENSG00000119408", "ENSG00000172057", "ENSG00000073605", "ENSG00000134262",
"ENSG00000149485", "ENSG00000134461", "ENSG00000176476", "ENSG00000134490",
"ENSG00000101384", "ENSG00000134242", "ENSG00000167740", "ENSG00000010030",
"ENSG00000159840", "ENSG00000213533", "ENSG00000118557", "ENSG00000112167",
"ENSG00000112139", "ENSG00000130338", "ENSG00000100321", "ENSG00000137757",
"ENSG00000146414", "ENSG00000135452", "ENSG00000140105", "ENSG00000100599",
"ENSG00000105656", "ENSG00000178952", "ENSG00000182108", "ENSG00000214113",
"ENSG00000111450", "ENSG00000073584", "ENSG00000055130", "ENSG00000103811",
"ENSG00000197728", "ENSG00000160226", "ENSG00000203760", "ENSG00000160194",
"ENSG00000196365", "ENSG00000166896", "ENSG00000013288", "ENSG00000168071",
"ENSG00000148824", "ENSG00000124532", "ENSG00000184293", "ENSG00000197982",
"ENSG00000167987", "ENSG00000175643", "ENSG00000184992", "ENSG00000132170",
"ENSG00000196850", "ENSG00000116793", "ENSG00000161395", "ENSG00000145029",
"ENSG00000180316", "ENSG00000167264", "ENSG00000198270", "ENSG00000075975",
"ENSG00000196449", "ENSG00000163453", "ENSG00000170191", "ENSG00000146192",
"ENSG00000065491", "ENSG00000183255", "ENSG00000197601", "ENSG00000198755",
"ENSG00000181004", "ENSG00000160469", "ENSG00000112599", "ENSG00000089327",
"ENSG00000051596", "ENSG00000161405", "ENSG00000104972", "ENSG00000130511",
"ENSG00000168685", "ENSG00000197119", "ENSG00000116977", "ENSG00000184076",
"ENSG00000066336", "ENSG00000100296", "ENSG00000188603", "ENSG00000035403",
"ENSG00000110665", "ENSG00000111231", "ENSG00000198040"), pchic = c("ENSG00000010810",
"ENSG00000013725", "ENSG00000057608", "ENSG00000061987", "ENSG00000065361",
"ENSG00000070985", "ENSG00000076321", "ENSG00000090104", "ENSG00000103811",
"ENSG00000105371", "ENSG00000105655", "ENSG00000111249", "ENSG00000111252",
"ENSG00000112182", "ENSG00000116741", "ENSG00000116747", "ENSG00000116750",
"ENSG00000116793", "ENSG00000117586", "ENSG00000117592", "ENSG00000117593",
"ENSG00000120334", "ENSG00000120337", "ENSG00000124356", "ENSG00000125245",
"ENSG00000125810", "ENSG00000126353", "ENSG00000130511", "ENSG00000134460",
"ENSG00000134882", "ENSG00000134954", "ENSG00000135439", "ENSG00000135452",
"ENSG00000135502", "ENSG00000136634", "ENSG00000140379", "ENSG00000141933",
"ENSG00000151655", "ENSG00000158246", "ENSG00000162892", "ENSG00000162894",
"ENSG00000162896", "ENSG00000162897", "ENSG00000166908", "ENSG00000168925",
"ENSG00000168928", "ENSG00000169508", "ENSG00000172243", "ENSG00000172322",
"ENSG00000172575", "ENSG00000173575", "ENSG00000175215", "ENSG00000177414",
"ENSG00000178498", "ENSG00000178726", "ENSG00000183734", "ENSG00000184281",
"ENSG00000185278", "ENSG00000185442", "ENSG00000185811", "ENSG00000196352",
"ENSG00000196498", "ENSG00000197406", "ENSG00000197992", "ENSG00000199179",
"ENSG00000202034", "ENSG00000207719", "ENSG00000208317", "ENSG00000211574",
"ENSG00000214548", "ENSG00000220981", "ENSG00000221949", "ENSG00000224713",
"ENSG00000226455", "ENSG00000230483", "ENSG00000231560", "ENSG00000234741",
"ENSG00000238475", "ENSG00000238754", "ENSG00000240338", "ENSG00000240771",
"ENSG00000251817", "ENSG00000251908", "ENSG00000252122", "ENSG00000252619",
"ENSG00000254588", "ENSG00000256803", "ENSG00000257354", "ENSG00000257568",
"ENSG00000258498", "ENSG00000258919", "ENSG00000258922", "ENSG00000259003",
"ENSG00000259054", "ENSG00000266933", "ENSG00000267607", "ENSG00000272888"
))

## load in gtf file
GTF<-fread('/scratch/ob219/DATA/JAVIERRE_GWAS/support/Homo_sapiens.GRCh37.75.genes.gtf')
GTF[,ensg:=gsub('.*"(.*)"',"\\1",tstrsplit(GTF$V9,';',fixed=TRUE,fill=TRUE)[[1]])]
GTF<-GTF[GTF$ensg %in% unlist(foo) & V2=='protein_coding',]

GTF$pchic<-GTF$ensg %in% foo$pchic
GTF$coloc<-GTF$ensg %in% foo$coloc
GTF$twas<-GTF$ensg %in% foo$twas



## bin ensg by biotype

## what is the distribution of non protein coding genes across the genome ?

library(rtracklayer)
t1d.gr<-import.bed(con='https://www.immunobase.org/downloads/regions-files-archives/latest_1.11/Hs_GRCh37-SLE-assoc_tableBED')
seqlevels(t1d.gr)<-gsub('chr','',seqlevels(t1d.gr))
seqlevels(t1d.gr)[seqlevels(t1d.gr)=='X']<-23


res<-lapply(c('coloc','pchic','twas'),function(f){
  message(f)
  max.coord<-GTF[,list(mc=max(V5)),by=V1]
  reso<-1e7
  genome<-lapply(split(max.coord,max.coord$V1),function(m){
      s<-seq(1,m$mc,by=reso)
      s2<-c(head(s+reso,-1)-1,m$mc)
      GRanges(seqnames=Rle(m$V1),ranges=IRanges(start=s,end=s2))
  })
  genome<-unlist(GRangesList(genome))
  GTF.pc.nc<-subset(GTF,V2=='protein_coding' & GTF[[f]])
  GTF.pc.nc[,coord:=V4]
  GTF.pc.nc[GTF.pc.nc$V7=='-',coord:=V5]
  neg<-GTF.pc.nc[GTF.pc.nc$V7=='-',]
  ps<-GTF.pc.nc[GTF.pc.nc$V7=='+',]
  gr<-rbind(ps[ps[, .I[which.min(V4)], by=ensg]$V1],neg[neg[, .I[which.min(V4)], by=ensg]$V1])[,.(V1,coord)]
  gr<-with(gr,GRanges(seqnames=Rle(V1),range=IRanges(start=coord,end=coord)))
  ol<-as.matrix(findOverlaps(genome,gr))
  counts<-sapply(split(ol[,2],ol[,1]),length)
  genome$gc<-0
  genome[as.numeric(names(counts)),]$gc<-counts
  ## next annotate genome regions that overlap susceptibility regions
  ol<-as.matrix(findOverlaps(genome,t1d.gr))
  genome$t1d<-FALSE
  if(nrow(ol)!=0)
    genome[ol[,1],]$t1d<-TRUE
  gt<-data.table(data.frame(genome))
  gt$chr<-as.numeric(as.character(gt$seqnames))
  gt[gt$seqnames=='X']$chr<-23
  gt[gt$seqnames=='Y']$chr<-24
  gt<-gt[order(gt$chr,gt$start),]

  tmp<-split(gt$end,gt$chr)
  tmp2<-split(gt$start,gt$chr)

  cs<-c(0,head(cumsum(as.numeric(sapply(tmp,max))),-1))+1
  for(i in seq_along(tmp)){
    tmp2[[i]]<-tmp2[[i]] + cs[i]
  }

  gt$pstart<-do.call('c',tmp2)
  gt[,pend:=pstart+width]
  gt$fac<-f
  gt
})


res<-rbindlist(res)
res<-res[res$gc!=0,]

library(ggplot2)

ggplot(res,aes(xmin=pstart,xmax=pend,ymin=0,ymax=gc,fill=t1d)) + geom_rect() + theme_bw() + facet_grid(fac~.)

ggplot(gt,aes(xmin=pstart,xmax=pend,ymin=0,ymax=gc,fill=as.factor(chr%%2))) + geom_rect() + theme_bw()
break()
ggplot(gt,aes(x=pstart,y=gc,color=as.factor(chr%%2))) + geom_path() + theme_bw()


## what is the average size of captured sequence

round(mean(width(with(DT,IRanges(start=V2,end=V3)))))

## what is the average size for all HindIII for the genome

h3<-fread('/scratch/ob219/DATA/JAVIERRE_GWAS/support/Digest_Human_HindIII.bed')
round(mean(width(with(h3,IRanges(start=V2,end=V3)))))


## protein coding promoters captured per gene
DT.pc<-DT[DT$bt=='protein_coding',]
by.frag<-split(DT.pc$ensg,DT.pc$V4)
c<-data.table(table(sapply(by.frag,function(f) length(unique(f)))))
(1-(c[1,]$N/length(by.frag))) * 100

library(magrittr)
library(data.table)
library(biomaRt)
library(pheatmap)
library(parallel)
library(ggplot2)

DATA.DIR<-'/Users/oliver/DATA/JAVIERRE_GWAS/out/geneScore'
fs<-list.files(path=DATA.DIR,pattern="*.tab",full.names = TRUE)
all.dat<-rbindlist(lapply(fs,fread))
all.dat.pc<-all.dat[all.dat$biotype=='protein_coding',]

gs <- all.dat.pc[,.(disease,ensg=ensg,gs=all_gene_score)]
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','gene_biotype','external_gene_name','entrezgene'), filters = 'ensembl_gene_id', values = unique(gs$ensg), mart =ensembl) %>% data.table
M<-merge(gs,genedesc,by.x='ensg',by.y='ensembl_gene_id',all.x=TRUE)[!is.na(entrezgene),]
dupent <-  M[which(duplicated(M[,.(disease,entrezgene)])),]$entrezgene %>% unique
M <- M[!entrezgene %in% dupent,]
M <- M[disease!='RA',]

tab<-fread("/Users/oliver/DATA/JAVIERRE_GWAS/support/gwas_manifest.csv")
M<-merge(M,tab,by.x='disease',by.y='label')
M[disease=='RA_OKADA_IMB',disease:='RA']
M[disease=='CD_IMB',disease:='CD']
M[,disease_lab:=gsub("\\_"," ",disease)]
M[disease=='T2D',category:='Other']
levs <- M[,list(cat=unique(category)),by=c('disease_lab','disease')][order(cat,disease),]$disease_lab
M[,disease_lab:=factor(disease_lab,levels=levs)]

library(cowplot)
library(ggridges)

ppd <- ggplot(M,aes(x=gs,y=disease_lab,fill=category)) + geom_density_ridges(col="grey") + coord_cartesian(xlim=c(0,0.05)) +
scale_fill_manual(name = "Trait Category",values = c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")) +
xlab("Overall COGS Score") + ylab("Trait/COGS score density")

save_plot("~/tmp/cogs_score_distro.pdf",ppd,base_width=8,base_height=7)

## what is coverage between traits ?
miss.genes <- M[,list(count=.N),by=entrezgene][count!=31,]$entrezgene
## to make the universe the same across all traits remove genes that are not covered by all traits

## load in hallmark
gmt <- scan("~/tmp/h.all.v6.1.entrez.gmt",character(),sep="\n")

## load in reactome
gmt <- scan("~/tmp/c2.cp.reactome.v6.2.entrez.gmt",character(),sep="\n")
hm <- list()
for(i in 1:length(gmt)){
  v <- strsplit(gmt[i],"\t") %>% unlist
  hm[[v[1]]] <- as.integer(v[3:length(v)])
}
names(hm) <- gsub("REACTOME_","",names(hm))
#names(hm) <- gsub("HALLMARK_","",names(hm))

## try a fisher.test

ft <- mclapply(unique(M$disease_lab) %>% as.character,function(d){
  message(sprintf("Processing %s",d))
  t.DT <- M[disease_lab==d & !is.na(gs),]
  gsv <- t.DT$gs>0.5
  names(gsv) <- t.DT$entrezgene
  ## create universe
  ## remove those that are not in universe of all genes
  adj.hm <- lapply(hm,function(x){
    x[x %in% names(gsv)]
  })
  pall <- lapply(adj.hm,function(gs){
    FT <- fisher.test(table(gsv,names(gsv) %in% gs))
    data.table(stat=FT$estimate,p=FT$p.value)
  }) %>% rbindlist
  pall[,c('disease','pw'):=list(d,names(adj.hm))]
  pall
},mc.cores=8) %>% rbindlist
ft[,p.adj:=p.adjust(p,method="fdr")]
ft[,p.adj.mlp:=-log10(p.adj)]
ft[,p.mlp:=-log10(p)]
mat <- melt(ft,id.vars=c('disease','pw'),measure.vars='p.adj.mlp') %>% dcast(.,disease~pw+variable)
rnames <- mat[[1]]
mat <- as.matrix(mat[,-1])
rownames(mat) <- rnames
colnames(mat) <- gsub("_p.adj.mlp","",colnames(mat))
## filter out those diseases with nothing going on !
sig.disease <- apply(mat,1,function(x) sum(x> -log10(0.05))>0)
sig.disease <- names(sig.disease[sig.disease])
sig.path <-  apply(mat,2,function(x) sum(x> -log10(0.05))>0)
sig.path <- names(sig.path[sig.path])
pheatmap(mat[sig.disease,sig.path])

## can we make this in ggplot so it looks nicer and we have more control. Clustering not that important apart from by disease

hcd <- dist(mat[sig.disease,sig.path]) %>% hclust
dorder <- sig.disease[hcd$order]
## actually order by category
dorder <- M[,list(cat=unique(category)),by=c('disease_lab','disease')][order(cat,disease),][disease_lab %in% sig.disease,]$disease_lab %>% as.character
pwcd <- t(mat[sig.disease,sig.path]) %>% dist %>% hclust
pworder <- sig.path[pwcd$order]
pwlabs <- gsub("\\_"," ",pworder)
pwlabs[grep("MEGA",pwlabs)] <- "MEGAKARYOCYTE DEVELOPMENT AND PLATELET PRODUCTION"
pwlabs[grep("ZAP",pwlabs)] <- "TRANSLOCATION OF ZAP 70 TO IMMUNOLOGICAL SYNAPSE"
#pworder <- sig.path
pft <- ft[disease %in% sig.disease & pw %in% sig.path,]
pft[,c('dfact','pwfact'):=list(factor(disease,levels=dorder),factor(pw,levels=pworder,labels=pwlabs))]
pft[,plot.mlp:=ifelse(p.adj<0.05,p.mlp,NA)]

cat.col <- c(Autoimmune="darkblue",Blood="red",Metabolic="orange",Other="black")

dorder.col <-  cat.col[M[,list(cat=unique(category)),by=c('disease_lab','disease')][order(cat,disease),][disease_lab %in% sig.disease,]$cat]

ftp <- ggplot(pft,aes(y=dfact,x=pwfact,fill=plot.mlp)) + geom_tile(col='black') +
scale_fill_continuous(name="-log10(P.adj)",low="thistle2", high="darkred", guide="colorbar",na.value="white") +
theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),axis.text.y = element_text(colour = dorder.col)) + xlab("Reactome Pathway") + ylab("Trait")

save_plot("~/tmp/fisher_test_cogs.pdf",ftp,base_width=12,base_height=8)

# dl <- mclapply(unique(M$disease),function(d){
#   message(sprintf("Processing %s",d))
#   t.DT <- M[disease==d & !is.na(gs),]
#   gsv <- t.DT$gs
#   names(gsv) <- t.DT$entrezgene
#   ## create universe
#   ## remove those that are not in universe of all genes
#   adj.hm <- lapply(hm,function(x){
#     x[x %in% names(gsv)]
#   })
#   pall <- lapply(adj.hm,function(gs){
#     inset <- which(names(gsv) %in% gs)
#     W <- wilcox.test(gsv[inset], gsv[-inset])
#     data.table(stat=W$statistic,p=W$p.value)
#   }) %>% rbindlist
#   pall[,c('disease','pw'):=list(d,names(adj.hm))]
#   pall
# },mc.cores=8) %>% rbindlist
#
# dl[,p.adj:=p.adjust(p,method="fdr")]
# dl[,p.adj.mlp:=-log10(p.adj)]
# dl[,p.mlp:=-log10(p)]
# sig.path<-unique(dl[p.adj<0.01,]$pw)
#
# mat <- melt(dl,id.vars=c('disease','pw'),measure.vars='p.adj.mlp') %>% dcast(.,disease~pw+variable)
# rnames <- mat[[1]]
# mat <- as.matrix(mat[,-1])
# rownames(mat) <- rnames
# colnames(mat) <- gsub("_p.adj.mlp","",colnames(mat))


#pheatmap(mat[,sig.path])

dl.t <- mclapply(unique(M$disease),function(d){
  message(sprintf("Processing %s",d))
  t.DT <- M[disease==d & !is.na(gs),]
  gsv <- t.DT$gs
  names(gsv) <- t.DT$entrezgene
  ## create universe
  ## remove those that are not in universe of all genes
  adj.hm <- lapply(hm,function(x){
    x[x %in% names(gsv)]
  })
  pall <- lapply(hm,function(gs){
      inset <- which(names(gsv) %in% gs)
      K <- ks.test(gsv[inset], gsv[-inset])
      data.table(stat=K$statistic,p=K$p.value)
    }) %>% rbindlist
    pall[,c('disease','pw'):=list(d,names(adj.hm))]
    pall
},mc.cores=8) %>% rbindlist

dl.t[,p.adj:=p.adjust(p,method="fdr")]
dl.t[,p.adj.mlp:=-log10(p.adj)]
dl.t[,p.mlp:=-log10(p)]
sig.path<-unique(dl.t[p.adj<0.01,]$pw)

mat <- melt(dl.t,id.vars=c('disease','pw'),measure.vars='p.adj.mlp') %>% dcast(.,disease~pw+variable)
rnames <- mat[[1]]
mat <- as.matrix(mat[,-1])
rownames(mat) <- rnames
colnames(mat) <- gsub("_p.adj.mlp","",colnames(mat))
pheatmap(mat[,sig.path])

## what about fisher.test

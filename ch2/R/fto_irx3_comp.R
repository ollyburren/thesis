gene2pubmed <- fread("~/tmp/gene2pubmed")

## we know that FTO is 79068
## IRX3 is 79191
## IRX5 is 10265

FTO.ids <- gene2pubmed[GeneID==79068,]$PubMed_ID
IRX3.ids <- gene2pubmed[GeneID==79191,]$PubMed_ID
IRX5.ids <- gene2pubmed[GeneID==10265,]$PubMed_ID

## use easy pubmed to get meta data about each publication

foo<-split(FTO.ids,1:length(FTO.ids))
names(foo) <- rep('Id',length(foo))

q<-list(IdList=foo)

library(rentrez)
library(parallel)


getDates <- function(pids){
  res <- lapply(split(FTO.ids, ceiling(seq_along(pids)/25)),function(ids){
    message("chunk")
    pxml<- entrez_fetch(db = "pubmed", id =ids,
                         rettype = "xml", parsed = T)
    foo<-lapply(XML::xmlToList(pxml),function(x) data.table(pmid=x$MedlineCitation$PMID$text,year=x$PubmedData$History$PubMedPubDate$Year))   %>% rbindlist
    return(foo)
  }) %>% rbindlist
  return(res)
}


search_year <- function(year, term){
    query <- paste(term, "AND (", year, "[PDAT])")
    entrez_search(db="pubmed", term=query, retmax=0)$count
}

year <- 2001:2017
papers.FTO <- sapply(year, search_year, term="FTO", USE.NAMES=FALSE)
papers.IRX3 <- sapply(year, search_year, term="IRX3", USE.NAMES=FALSE)
papers.IRX5 <- sapply(year, search_year, term="IRX5", USE.NAMES=FALSE)

FTO.DT <- data.table(years=year,papers=papers.FTO,Gene='FTO')
IRX3.DT <- data.table(years=year,papers=papers.IRX3,Gene='IRX3')
IRX5.DT <- data.table(years=year,papers=papers.IRX5,Gene='IRX5')

plot.DT <- list(FTO.DT,IRX3.DT,IRX5.DT) %>% rbindlist

library(cowplot)

pp <- ggplot(plot.DT[years<2016 & Gene !='IRX5',],aes(x=years,y=papers,color=Gene,group=Gene)) +
geom_point() + geom_line() + geom_vline(xintercept=2007,lty=2) + geom_vline(xintercept=2014,lty=2) +
xlab("Publication Year") + ylab("Publication Count")
save_plot("~/tmp/fto_pub_comparison.pdf",pp)

entrez_search(db="pubmed", term="FTO", retmax=0)$count

FTO.DT <- getDates(FTO.ids)
IRX3.DT <- getDates(IRX3.ids)
IRX5.DT <- getDates(IRX5.ids)

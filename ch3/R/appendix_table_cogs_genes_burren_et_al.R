## code to generate tables of different results for appendix
library(data.table)
library(magrittr)
DT <- fread("~/Downloads/13059_2017_1285_MOESM9_ESM.tab")
DT[analysis=='GWAS_SS',analysis:='GW']
DT[analysis=='ICHIP_GF',analysis:='IC']
unique(DT$Ensembl_GeneID) %>% length

## change the colnames

setnames(DT,c('ENSG','Name','Disease','Analysis','COGS','Context','COGS_DROP','High_DROP','Expression','eRNA','Sus_DROP','Locus','SNP','p-value'))

DT <- DT[,.(ENSG,Name,Disease,Analysis,COGS,Context,Expression,eRNA,Locus,SNP,`p-value`)]
DT[,COGS:=signif(COGS,digits=2)]
## recode Context
DT[Context=='interaction',Context:='INT']
DT[Context=='promoter',Context:='PROM']
DT[Context=='Total_CD4_Activated',Context:='ACT']
DT[Context=='Total_CD4_NonActivated',Context:='NACT']
DT[Context=='coding',Context:='CODE']
DT[Context=='noncoding',Context:='NCODE']
DT[Context=='overall',Context:='OVERALL']

DT[is.na(Expression),Expression:='ND']
DT[Expression=='up',Expression:='+']
DT[Expression=='down',Expression:='-']
DT[Expression=='nsig',Expression:='=']

DT[is.na(eRNA),eRNA:='ND']
DT[eRNA=='up',eRNA:='+']
DT[eRNA=='down',eRNA:='-']
DT[eRNA=='nsig',eRNA:='=']

M[Locus=='NONE',Locus:='']
M[is.na(Locus),Locus:='']
M[SNP=='NONE',SNP:='']
M[is.na(SNP),SNP:='']



## get gene positions from biomart so we can sort these correctly.

library(biomaRt)
ensembl = useEnsembl(host="grch37.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
genedesc <- getBM(attributes=c('ensembl_gene_id','chromosome_name','start_position'),
filters = 'ensembl_gene_id', values = DT$ENSG, mart =ensembl) %>% data.table

M <- merge(DT,genedesc,by.x='ENSG',by.y='ensembl_gene_id')
M <- M[order(as.numeric(chromosome_name),start_position),.(ENSG,Name,Disease,Analysis,COGS,Context,Expr=Expression,eRNA,Locus,SNP,`p-value`)]
library(xtable)
digi <- rep(2,ncol(M)+1)
digi[12] <- -1
xtable(M,digits=digi)

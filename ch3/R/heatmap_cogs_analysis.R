## summary diagram using a heatmap across all nodes and all traits


## load in the data

#pf <- list.files(path="/home/ob219/rds/hpc-work/thesis/hgenescore_0.5",pattern='*.pmi_prioritised.tab',full.names=TRUE)
pf <- list.files(path="/home/ob219/rds/hpc-work/thesis/hgenescore_0.01",pattern='*.pmi_prioritised.tab',full.names=TRUE)


DT <- lapply(pf,fread) %>% rbindlist

DT <- DT[!disease %in% c('RA','T1D'),]
DT[disease=='RA_OKADA_IMB',disease:='RA']
DT[disease=='COOPER_T1D',disease:='T1D']
DT[disease=='CD_IMB',disease:='CRO']

DT[node=='Neutrophils',node:='Neu']
DT[node=='Endothelial_precursors',node:='EndP']
DT[node=='Monocytes',node:='Mon']
DT[node=='Total_CD4_MF',node:='tCD4']
DT[node=='Erythroblasts',node:='Ery']
DT[node=='Macrophages_M0',node:='M0']
DT[node=='Macrophages_M1',node:='M1']
DT[node=='Macrophages_M2',node:='M2']
DT[node=='Naive_CD4',node:='nCD4']
DT[node=='Foetal_thymus',node:='FetT']
DT[node=='Megakaryocytes',node:='MK']
DT[node=='Total_CD4_NonActivated',node:='naCD4']
DT[node=='Total_CD4_Activated',node:='aCD4']
DT[node=='Total_B',node:='tB']
DT[node=='Naive_B',node:='nB']
DT[node=='Naive_CD8',node:='nCD8']
DT[node=='Total_CD8',node:='tCD8']

DT[,category:='Other']
DT[disease %in% c('MS','CEL','T1D','CRO','PBC','UC','SLE','RA'),category:='Immune']
DT[disease %in% c('BMI','LDL','TG','HDL','TC','INS','INS_BMI','GLUCOSE','GLUCOSE_BMI'),category:='Metabolic']
DT[disease %in% c('HB','MCH','PCV','MCHC','RBC','MCV','PLT','PV'),category:='Blood']

disease.cat <- DT[,.(disease,category)] %>% unique %>% as.data.frame
rownames(disease.cat) <- disease.cat$disease
disease.cat$disease <- NULL
colnames(disease.cat) <- 'Trait Category'

tissue.cat <- DT[,.(node,isLeaf=ifelse(isLeaf==TRUE,'True','False'))] %>% unique %>% as.data.frame
rownames(tissue.cat) <- tissue.cat$node
tissue.cat$node <- NULL
colnames(tissue.cat) <- 'Is Leaf Node'

an.col <- list('Trait Category' = c('Blood'='red','Immune'='darkblue','Metabolic'='darkgreen','Other'='black'))



## remove non tissue specific nodes

#f.DT <- DT[!node %in% c('overall','noncoding','interaction','promoter','coding'),]
f.DT <- DT[!node %in% c('overall','noncoding','promoter','coding'),]

M <- melt(f.DT,id.vars='disease',measure.vars='node') %>% dcast(.,disease~variable+value)

## recode the leafnodes so they are consistent with table 2.1




M.mat <- as.matrix(M[,-1])
rownames(M.mat) <- M$disease
cnames <- colnames(M.mat) %>% gsub("node\\_","",.)
cnames[cnames=='Total_CD4_MF'] <- 'Total_CD4'
colnames(M.mat) <- cnames
library(pheatmap)
#pheatmap(log10(M.mat + 1),annotation_row=disease.cat,display_numbers=M.mat,annotation_col=tissue.cat,annotation_colors=an.col,clustering_method='ward.D2',filename="~/tmp/hcogs_results.pdf",legend=FALSE,width=8)
pheatmap(log10(M.mat + 1),annotation_row=disease.cat,annotation_col=tissue.cat,annotation_colors=an.col,clustering_method='ward.D2',filename="~/tmp/hcogs_results_0.01.pdf",legend=FALSE,width=8)


## get some stats what is the IQR of counts assigned to leaf nodes

M.mat[,colnames(M.mat) %in% rownames(tissue.cat)[tissue.cat[[1]]=='True']] %>% as.numeric %>% summary
M.mat[,!colnames(M.mat) %in% rownames(tissue.cat)[tissue.cat[[1]]=='True']] %>% as.numeric %>% summary

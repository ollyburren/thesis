## modeling sharing between beta's for thesis

bcor <- function(scase,sctrl,n11,n12,n01,n02){
  grat <- sqrt((n11 * n12) / (n01 * n02))
  gratv <- sqrt((n01 * n02) / (n11 * n12))
  tot <- (sctrl * grat) + (scase * gratv)
  sprintf("grat %f, gratv %f, tot %f",grat,gratv,tot) %>% print
  tot/sqrt((n11 + n01) * (n12 + n02))
}

## for crohn's disease
bcase <- 12294
bctrl <- 28072
ucase <- 1032
uctrl <- 336127

case_share <- seq(0,1,by=0.1)
ctrl_share <- seq(0,1,by=0.1)
## this is the sharing for crohn's disease
bcor(min(ucase,bcase) * case_share, min(bctrl,uctrl) * 0,bcase,ucase,bctrl,uctrl)

bcor(min(ucase,bcase) * seq(0,1,by=0.1), min(bctrl,uctrl) * 1,bcase,ucase,bctrl,uctrl)

bcor(min(ucase,bcase) * case_share,min(bctrl,uctrl) * 0.5,bcase,ucase,bctrl,uctrl)

dodat <- lapply(case_share,function(cs){
  #lapply(ctrl_share,function(ctrs){
    bc <- bcor(min(ucase,bcase) * cs,min(bctrl,uctrl) * ctrl_share,bcase,ucase,bctrl,uctrl)
    data.table(case_share=cs,ctrl_share=ctrl_share,beta_cor=bc)
  #}) %>% rbindlist
}) %>% rbindlist

dodat[,ctrl_share:=factor(ctrl_share)]
ggplot(dodat,aes(x=case_share,y=beta_cor,col=ctrl_share)) + geom_point()


dodat <- lapply(case_share,function(cs){
  #lapply(ctrl_share,function(ctrs){
    bc <- bcor(min(bcase,bcase) * cs,min(bctrl,bctrl) * ctrl_share,bcase,bcase,bctrl,bctrl)
    data.table(case_share=cs,ctrl_share=ctrl_share,beta_cor=bc)
  #}) %>% rbindlist
}) %>% rbindlist

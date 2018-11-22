library(ggplot2)
library(data.table)
library(cowplot)

#sets the seed for random number generation.
set.seed(41)

n<-100
x <- 1:n
#100 normally distributed rand. nos.
ex <- rnorm(n, 0, 10^2)
ey <- rnorm(n, 0, 10^2)
# y is a linear function of x
y <- 30 + 10 * x
#adds "noise" to x and y
x_obs <- x + ex
y_obs <- y + ey
P <- data.table(x=x_obs,y=y_obs)
#compute the centroid of P
centroid<-data.table(x=mean(P$x),y=mean(P$y))
# next centre
M <- cbind(x_obs-centroid$x,y_obs-centroid$y)
## covariance matrix is
Mcov <- cov(M)
## alternative way of doing this is for linear algebra understanding
SM <- scale(M,center = TRUE,scale=FALSE) # get deviance (i.e. subtract mean of the cols)
(t(SM) %*% SM)/(nrow(SM)-1)

## get eigen vectors for covariance matrix
eigenVal <- eigen(Mcov)$value
eigenVectors <- eigen(Mcov)$vectors
# we can compute first principal axis as follows
pc1 <- data.table(x=x_obs,y=eigenVectors[2,1]/eigenVectors[1,1]*M[x]+mean(y_obs))
# this is the second
pc2 <- data.table(x=x_obs,y=eigenVectors[2,2]/eigenVectors[1,2]*M[x]+mean(y_obs))

## PLOT 2 show projection of data into PC1
trans <- (M%*%eigenVectors[,1])%*%eigenVectors[,1] #compute projections of points
P_proj <- scale(trans, center=-cbind(mean(x_obs),mean(y_obs)), scale=FALSE)
P_proj <- data.table(x=P_proj[,1],y=P_proj[,2])
plines <- cbind(x_obs,y_obs,P_proj)
setnames(plines,c('x','y','p.x','p.y'))

## for this projection points 18 and 99 are the most extreme for each PC so we want
# to hilight
labels<-character(length=100)
labels[c(18,99)] <- c(18,99)

hilight <- rep('dodgerblue',100)
hilight[c(18,99)] <- 'red'

size <- rep(0.5,100)
size[c(18,99)] <- 2


# to keep square makes sure x and y axis are the same
limsits <- c(min(c(x_obs,y_obs)),max(c(x_obs,y_obs)))

# P

ppa<-ggplot(P,aes(x=x,y=y)) + geom_point() +
xlab("X") + ylab("Y") + coord_cartesian(xlim=limsits,ylim=limsits) +
geom_text(label=labels,hjust=-0.5) + geom_vline(xintercept=centroid$x,alpha=0.1) +
geom_hline(yintercept=centroid$y,alpha=0.1) + geom_point(data=centroid,col='firebrick',size=3) +
geom_line(data=data.table(x=limsits,y=limsits+30),alpha=0.5,lty=2)

## get a line passing through centoid with a random intercept


library(magrittr)


projPoint <- function(x,y,a,b,c){
  x1 = ((b * (b * x - a * y)) - (a * c))/(a^2+b^2)
  y1 = ((a * (-1 * b * x + a * y)) - (b * c))/(a^2+b^2)
  data.table(x=x1,y=y1)
}

simLine <- function(x,y,gradient){
  mx <- mean(x)
  my <- mean(y)
  intercept <- my + -1 * gradient * mx
  b<- -1
  proj <- projPoint(x,y,gradient,b,intercept)
  #dist <- sum(sqrt((proj$x-mx)^2 + (proj$y-my)^2))
  dist <- sum(sqrt((proj$x-x)^2 + (proj$y-y)^2))
  var <- mean(sqrt((proj$x-mx)^2 + (proj$y-my)^2))
  #var <- var(proj$x) + var(proj$y) + 2 * cov(proj$x,proj$y)
  #var <- cov(proj$x,proj$y)
  data.table(sum.distance=dist,variance=var,grad = gradient)
}

projP <- function(x,y,gradient){
  mx <- mean(x)
  my <- mean(y)
  intercept <- my + -1 * gradient * mx
  b<- -1
  projPoint(x,y,gradient,b,intercept)
}

foo <- cbind(P,projP(P$x,P$y,1))
setnames(foo,c('x','y','x1','y1'))

#foo<-rbind(P[,type:="raw"],projP(P$x,P$y,100)[,type:='proj'])

ggplot(P,aes(x=x,y=y)) + geom_point() + geom_smooth(method="lm",formula="y~x") + geom_point(data=centroid,col='red',size=3)
#geom_segment(data=foo,aes(xend=x1,yend=y1),col='red',lty=2,alpha=0.5) + geom_smooth(method="lm",formula="y~x")

ggplot(res,aes(x=sum.distance,y=variance)) + geom_point()

res <- lapply(seq(0,10,by=0.01),function(grad){
  simLine(P$x,P$y,grad)
}) %>% rbindlist

ggplot(res,aes(x=sum.distance,y=variance)) + geom_point()
ggplot(res,aes(x=grad,y=sum.distance)) + geom_point()
ggplot(res,aes(x=grad,y=variance)) + geom_point()
## compute the distance and variance


ppa<-ggplot(P,aes(x=x,y=y)) + geom_point() + geom_smooth(method="lm") +
xlab("X") + ylab("Y") + coord_cartesian(xlim=limsits,ylim=limsits) +
geom_text(label=labels,hjust=-0.5) + geom_vline(xintercept=centroid$x,alpha=0.1) +
geom_hline(yintercept=centroid$y,alpha=0.1) + geom_point(data=centroid,col='firebrick',size=3) +
geom_line(data=data.table(x=limsits,y=limsits+30),alpha=0.5,lty=2) + geom_point(data=lineSim(P$x,-82),col='green')



ppb<-ggplot(P,aes(x=x,y=y)) + geom_point(col='grey') +
xlab("X") + ylab("Y") + geom_line(data=pc1,col='black') +
geom_line(data=pc2,col='grey') + coord_cartesian(xlim=limsits,ylim=limsits) + geom_point(data=P_proj,col=hilight,size=size) +
geom_segment(data=plines,aes(xend=p.x,yend=p.y),col=hilight,lty=2,alpha=0.5) + geom_point(data=centroid,col='firebrick',size=3) +
geom_text(label=labels,hjust=-0.5) + geom_point(data=lineSim(P$x,-82),col='green')

trans <- (M%*%eigenVectors[,2])%*%eigenVectors[,2] #compute projections of points
P_proj2 <- scale(trans, center=-cbind(mean(x_obs),mean(y_obs)), scale=FALSE)
P_proj2 <- data.table(x=P_proj2[,1],y=P_proj2[,2])
plines2 <- cbind(x_obs,y_obs,P_proj2)
setnames(plines2,c('x','y','p.x','p.y'))
plines2[,l.x:=sqrt((x-p.x-mean(x_obs))^2 + (y-p.y-mean(y_obs))^2)]

ppc<-ggplot(P,aes(x=x,y=y)) + geom_point(col='grey') +
xlab("X") + ylab("Y") + geom_line(data=pc1,col='grey') +
geom_line(data=pc2,col='black') + coord_cartesian(xlim=limsits,ylim=limsits) + geom_point(data=P_proj2,,col=hilight,size=size) +
geom_segment(data=plines2,aes(xend=p.x,yend=p.y),col=hilight,lty=2,alpha=0.5) + geom_point(data=centroid,col='firebrick',size=3) +
geom_text(label=labels,hjust=-0.5)


## loadings in the state can be computed through pythagoras theorem
P_proj[,p.l:=sqrt((x-mean(x_obs))^2 + (x-mean(x_obs))^2) * sign(x-mean(x_obs))]
P_proj2[,p.l:=sqrt((x-mean(x_obs))^2 + (x-mean(x_obs))^2) * sign(x-mean(x_obs))]

ppd <- ggplot(data.table(x=P_proj$p.l,y=P_proj2$p.l),aes(x=x,y=y)) + geom_point() +
geom_text(label=labels,vjust=-1) + xlab("PC1 Loading") + ylab("PC2 Loading") +
geom_vline(xintercept=0,alpha=0.1) +
geom_hline(yintercept=0,alpha=0.1) + geom_point(data=data.table(x=0,y=0),col='firebrick',size=3)


p<-plot_grid(ppa,ppb,ppc,ppd,labels=c('A','B','C','D'))
save_plot("~/tmp/pca.pdf",p,base_width=10,base_height=10)

## what happens to distance between points and variance of points




# ## example PCA plot for region strat - do this better using 1K genomes data for different ethnicities.
# library(snpStats)
# library(data.table)
# library(ggplot2)
# snps <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1900_1900_500/snps.RDS')
# samples <- readRDS('/home/ob219/rds/hpc-work/as_basis/t1d_bootstrap/sims/wtccc_1900_1900_500/samples.RDS')
#
# cont.idx <- which(samples$cc==0)
#
# csnps <- snps[cont.idx,]
#
# cs <- apply(csnps,2,as.numeric)
# csf<-cs[,which(apply(cs,2,function(x) sum(x>3)) == 0)]
#
# pc<-prcomp(csf,scale=TRUE)
#
# ## downsample to 1000 snps
#
# pc.load <- cbind(data.table(pc$x),samples[cont.idx,]$b58region)
# ggplot(pc.load,aes(x=PC1,y=PC2,col=as.factor(V2))) + geom_point()

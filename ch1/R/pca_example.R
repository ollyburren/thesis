library(ggplot2)
library(data.table)
library(cowplot)

set.seed(42)             #sets the seed for random number generation.

n<-100
x <- 1:n              #creates a vector x with numbers from 1 to 100
ex <- rnorm(n, 0, 30) #100 normally distributed rand. nos. w/ mean=0, s.d.=30
ey <- rnorm(n, 0, 30) # " "
#y <- 30 + 2 * x
y <- 30 +  x         #sets y to be a vector that is a linear function of x
x_obs <- x + ex         #adds "noise" to x
y_obs <- y + ey         #adds "noise" to y
P <- data.table(x=x_obs,y=y_obs) #places points in matrix
centroid<-data.table(x=mean(P$x),y=mean(P$y))
# next centre
M <- cbind(x_obs-centroid$x,y_obs-centroid$y)
Mcov <- cov(M)

eigenVal <- eigen(Mcov)$value
eigenVectors <- eigen(Mcov)$vectors
# we can compute first principal axis as follows
pc1 <- data.table(x=x_obs,y=eigenVectors[2,1]/eigenVectors[1,1]*M[x]+mean(y_obs))
pc2 <- data.table(x=x_obs,y=eigenVectors[2,2]/eigenVectors[1,2]*M[x]+mean(y_obs))
d <- svd(M)$d
v <- svd(M)$v
trans <- (M%*%v[,1])%*%v[,1] #compute projections of points
P_proj <- scale(trans, center=-cbind(mean(x_obs),mean(y_obs)), scale=FALSE)
P_proj <- data.table(x=P_proj[,1],y=P_proj[,2])
plines <- cbind(x_obs,y_obs,P_proj)
setnames(plines,c('x','y','p.x','p.y'))
# to keep square makes sure x and y axis are the same
limsits <- c(min(c(x_obs,y_obs)),max(c(x_obs,y_obs)))

ppa<-ggplot(P,aes(x=x,y=y)) + geom_point(col='grey') +
xlab("X") + ylab("Y") + geom_line(data=pc1,col='black') +
geom_line(data=pc2,col='grey') + coord_cartesian(xlim=limsits,ylim=limsits) + geom_point(data=P_proj,col='dodgerblue',size=0.7) +
geom_segment(data=plines,aes(xend=p.x,yend=p.y),col='dodgerblue',lty=2,alpha=0.5) + geom_point(data=centroid,col='firebrick',size=3)

trans <- (M%*%v[,2])%*%v[,2] #compute projections of points
P_proj <- scale(trans, center=-cbind(mean(x_obs),mean(y_obs)), scale=FALSE)
P_proj <- data.table(x=P_proj[,1],y=P_proj[,2])
plines <- cbind(x_obs,y_obs,P_proj)
setnames(plines,c('x','y','p.x','p.y'))

ppb<-ggplot(P,aes(x=x,y=y)) + geom_point(col='grey') +
xlab("X") + ylab("Y") + geom_line(data=pc1,col='grey') +
geom_line(data=pc2,col='black') + coord_cartesian(xlim=limsits,ylim=limsits) + geom_point(data=P_proj,col='dodgerblue',size=0.7) +
geom_segment(data=plines,aes(xend=p.x,yend=p.y),col='dodgerblue',lty=2,alpha=0.5) + geom_point(data=centroid,col='firebrick',size=3)

p<-plot_grid(ppa,ppb,labels=c('A','B'))
save_plot("~/tmp/pca.pdf",p,base_width=8)

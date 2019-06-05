read.table('resultados.dat')-> d
ra <- d$V1
dec<-d$V2
z<-d$V3
Mr <- d$V6
UR <- d$V15
GR <- d$V16
C <- d$V14
fdev <- d$V19
r50 <- d$V17
mu <- d$V21

z2<-subset(z, Mr < -18.05)
Mr2 <- subset(Mr, Mr< -18.05)
UR2 <- subset(UR, Mr< -18.05)
GR2 <- subset(GR, Mr< -18.05)
C2 <- subset(C, Mr< -18.05)
fdev2 <- subset(fdev, Mr< -18.05)
r502 <- subset(r50, Mr< -18.05)
mu2 <- subset(mu, Mr< -18.05)

u_r <- hist(UR2,breaks=seq(-200,200,.1),plot=FALSE)
g_r <- hist(GR2,breaks=seq(-200,200,.1),plot=FALSE)

#library(repr)
library(astro)
#options(repr.plot.width=10, repr.plot.height=5)
png('hist_color2.png')
par(family='serif')
par(lwd=1)
par(cex=1)
par(mar=c(5,5,1,1))
par(mfrow=c(1,2))
par(mgp=c(1.2,.4,0))
aplot(u_r$mids,u_r$count/sum(u_r$counts),type='s',xlim=c(1,3.4),xlab='u-r',ylab='',col='blue',lwd=2.2)
aplot(g_r$mids,g_r$count/sum(g_r$counts),type='s',xlim=c(0,1.3),xlab='g-r',ylab='',col='red',lwd=2.2)
dev.off()


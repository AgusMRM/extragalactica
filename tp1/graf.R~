read.table('resultados.dat')->d
ra <- d$V1
dec<-d$V2
z<-d$V3
Mr <- d$V6
UR <- d$V15
GR <- d$V16

library(repr)
options(repr.plot.width=4, repr.plot.height=2)
#attach(mtcars)
par(mfrow=c(1,2))
par(mar=c(5,2,1,1))
urh <- hist(UR,breaks=seq(-400,400,.1),freq=FALSE)
grh <- hist(GR,breaks=seq(-400,400,.1),freq=TRUE)
plot(urh$mids,urh$count,type='s',xlim=c(0.,4.),col='blue',xlab='U-B',ylab='',lwd=3)
plot(grh$mids,grh$count,type='s',xlim=c(-0.2,1.3),col='red',xlab='G-R',ylab='',lwd=3)

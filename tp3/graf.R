png('Kz1.png')
read.table('mgs.dat')->d
#alfa, delta, z, 5 mag petro, 5 enrojecimientos, 5 correciones k a .0, 5 correciones a .1
z <- d$V3
ku <-d$V19
kg <- d$V20
kr <- d$V21
ki <- d$V22
kz <- d$V23


par(oma=c(0.5,0.3,0.3,0.3))  
par(mar=c(5,5,1,2))     # margenes para el box,bottom,left,top,right
par(mgp=c(3,1,0))      # margenes for the axis title, axis,labels and axis line.
par(mgp=c(2,1,0)) 
par(mfrow=c(5,1))
plot(z,ku,pch=19,cex=.1,ylab=expression(K[u]))
plot(z,kg,pch=19,cex=.1,ylab=expression(K[g]))
plot(z,kr,pch=19,cex=.1,ylab=expression(K[r]))
plot(z,ki,pch=19,cex=.1,ylab=expression(K[i]))
plot(z,kz,pch=19,cex=.1,ylab=expression(K[z]))
dev.off()
mr <- d$V6-d$V11
#plot(z,mr,pch='.')

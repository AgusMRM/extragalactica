attach(mtcars)
par(mfrow=c(2,1))
par(oma=c(0.5,0.3,0.3,0.3))  
par(mar=c(5,5,1,2))     # margenes para el box,bottom,left,top,right
par(mgp=c(3,1,0))      # margenes for the axis title, axis,labels and axis line.
par(mgp=c(2,1,0))      
par(cex.axis=1,cex.lab=1)
par(family="serif")
par(lwd=1)
par(cex=1)

xrads <- seq(0, 2*pi, length=30)
sinx  <- sin(xrads)
cosx  <- cos(xrads)

plot(xrads,sinx,xlab='x',ylab='Sen(x)',type='b',col='orange',lwd=2,pch=19)

plot(xrads,cosx,xlab='x',ylab='Cos(x)',type='b',col='forestgreen',lwd=2,pch=15)

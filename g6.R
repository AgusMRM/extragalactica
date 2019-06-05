x <- rnorm(1000)  #distribucion normal random

hx <- hist(x,breaks=seq(-4,4,0.5))

#escribir hx en pantalla y ver
plot(hx$mids,hx$count,type='s',xlim=c(-3.5,3.5),main="",xlab="X",ylab="N(X)",,col="blue",lwd=2)

#Normalizado
#plot(hx$mids,hx$count/sum(hx$counts),type='s',xlim=c(-3.5,3.5),main="",xlab="X",ylab="N(X)",,col="blue",lwd=2,)

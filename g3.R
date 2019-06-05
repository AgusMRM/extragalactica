par(lwd=2)
par(mar=c(5,5,2,2))             # margenes exteriores
par(mgp=c(3.7,1.3,0))           #margenes interiores
par(cex.axis=1.2,cex.lab=1.3)   #expand ejes y labels
par(family="serif")
par(cex=1)                      #expand texto y simbolos

xrads <- seq(0, 2*pi, length=30)  #genero 30 puntos
sinx  <- sin(xrads)
#pch ptype

#Opcion 1: solo puntos
plot(xrads, sinx,col="royalblue",pch=19,xlab="x [Radians]",ylab="Sen(x)",type="p")

#Opcion 1: solo linea
#plot(xrads, sinx,col="royalblue",pch=19,xlab="x [Radians]",ylab="Sen(x)",type="l")

#Opcion 1: puntos y lineas
#plot(xrads, sinx,col="royalblue",pch=19,xlab="x [Radians]",ylab="Sen(x)",type="b")


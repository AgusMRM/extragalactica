read.table("tabla.dat",head=F) -> D

x <- D$V1
y <- D$V2
z <- D$V3

#Caso 1
plot(x,y)

#Caso2
#plot(x,y,xlim=c(0,6),xlab="x",ylab="Function",cex=2)

#Caso3
#plot(x,y,xlim=c(0,6),xlab="x",ylab="Function",cex=1,pch=15,col=rgb(1,0.5,0,0.5))

#Caso4
#plot(x,y,xlim=c(0,6),xlab="x",ylab="Function",col=rgb(1,0.5,0,1),type='l',lwd=3)

#Caso5
#plot(x,y,xlim=c(0,6),ylim=c(0,26),xlab="x",ylab="Function",cex=2,pch=19)
#points(x,z,cex=1,pch=20,col='red',type='b')

#Caso6
#plot(x,y,xlim=c(0,6),ylim=c(0,26),xlab="x",ylab="Function",cex=2,pch=19,col='turquoise4')
#points(x,z,cex=1,pch=21,col='red',type='b')
#legend(1,22,c("y","z"),pch=c(19,21),col=c("turquoise4","red"),lty=c(0,1))

#ajustar recta
lm(y~x) -> z
abline(z)
print(z)
#summary(z)   si escribo en la terminal aparecen los datos de la recta, mas detallados

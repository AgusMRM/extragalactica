read.table('datosgrafico.at')->d
plot(d$V4,d$V1,col='red',type='l',ylim=c(0,1),xlab=expression(paste(log,'(',rho,')')),ylab='fraccion de galaxias')
points(d$V4,d$V2,col='green',type='l',lty=2)
points(d$V4,d$V3,col='blue',type='l',lty=5)
legend(0,.8,lty=c(1,2,5),col=c('red','green','blue'),c('E','S0','S'))


#png('FLgalaxias.png')

read.table('hist.dat')->d

M <- d$V2
L <- d$V1

plot(M,log10(d$V1),type='p',lwd=2,xlab=expression(M[r]),ylab=expression(paste(log,'(', Phi, '(', M[r], ')', ')')   ))
fit <- nls(L~I(a1-0.4*(M-a2)*(a3+1)-10**(-0.4*(M-a2))/log(10)),start=list(a1=4.32,a2=-19.88,a3=-0.75))

lines(M,predict(fit),col='green',lty=1,type='l')

#dev.off()

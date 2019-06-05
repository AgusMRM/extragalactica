library(magicaxis)

xrads <- seq(0, 2*pi, length=30)
sinx  <- sin(xrads)

plot(x,y,axes=FALSE,xlab=expression(alpha))
magaxis(side=1, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE)
magaxis(side=2, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=TRUE,las=2)
magaxis(side=3:4, majorn=5, minorn=5, tcl=0.5, ratio=0.5, labels=FALSE)

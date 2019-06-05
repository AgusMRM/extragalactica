x <- rnorm(1000)  #distribucion normal random

hist(x,breaks=10, col="blue", main="",freq=F,density=40,border="cyan")	##freq=F histograma normalizado a Area=1
box()

#hist(x,breaks=seq(-4,4,0.5), col="blue", main="",freq=F,density=30,border="royalblue")

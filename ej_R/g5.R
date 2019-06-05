
x<- c(1:5)
y<- c(1.1, 1.5, 2.9, 3.8, 5.2)
sd<- c(0.2, 0.3, 0.2, 0.1, 0.4)


plot (x, y, ylim=c(0, 6),pch=19)
#polygon(c(x,rev(x)),c(y+sd,rev(y-sd)),col=rgb(0.8,0.8,0.8,0.5),border=NA)
segments(x, y-sd,x, y+sd)

epsilon = 0.02
segments(x-epsilon,y-sd,x+epsilon,y-sd)
segments(x-epsilon,y+sd,x+epsilon,y+sd)



xx <- data.frame(x,y,sd)
print(xx)

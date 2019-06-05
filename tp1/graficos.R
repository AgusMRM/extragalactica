read.table('tabla2.csv',sep=',') -> d
z <- d$V3
mr <- d$V7
er <- d$V17

plot(z,(mr-er),pch='.')

x <- subset(z,mr-er>=14.5 & mr-er<17.77)
y <- subset(mr-er,mr-er>=14.5 & mr-er<17.77)
points(x,y,pch='.',col='red')

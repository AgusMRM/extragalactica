read.table('tabla2.csv',sep=',')-> D
z <- D$V3
hist(z,breaks=seq(0.02,.05,.0005))



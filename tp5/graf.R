read.table('field.dat') -> D

#Mr <- subset(D$V4, D$V21 > 0 & D$V21 < 20)
#nuvr <- subset(D$V21, D$V21 > 0 & D$V21 < 20)

#plot(Mr,nuvr,pch='.',col='forestgreen',ylim=c(0,8),xlim=c(-23,-16))
#abline(h=4)
#abline(h=5)

dn <- subset(D$V12, D$V12 < 4  & D$V12 > 1.6 & D$V9 > 8 & D$V12 >-90)
#ur <- subset(D$V5, D$V12 < 4 & D$V12 > 0  & D$V9 > 8 & D$V12 >-90)

#plot(dn,ur,pch='.',ylim=c(0,4),xlim=c(.5,3))
 

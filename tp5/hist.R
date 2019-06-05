read.table('field.dat') ->D
dn <- subset(D$V12, D$V12 > 0 & D$V12 < 3)

hist(dn,breaks=20)
abline(v=1.6)

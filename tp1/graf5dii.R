read.table('resultados.dat')->d
Mr <- d$V6
R <- d$V17
C <- d$V14

mr <- subset(Mr, c < 2.55)
r <- subset(R, c < 2.55)

plot(log10(r),Mr)

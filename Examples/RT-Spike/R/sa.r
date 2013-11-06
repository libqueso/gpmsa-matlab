me.CT <- as.matrix(read.table("../ME.CT"))
me.CB <- as.matrix(read.table("../ME.CB"))
me.CD <- as.matrix(read.table("../ME.CD"))

yl <- range(c(as.vector(me.CT),as.vector(me.CB),as.vector(me.CD)))
postscript("spike_sa.ps",height=3,width=8,paper="special",horizontal=F)
par(mfrow=c(1,3),pty="s",mgp=c(2,1,0),oma=c(0,0,0,0),
    mar=c(3,3,1,1))
palette(rainbow(ncol(me.CT)))
plot(0,0,xlab="Atwood number",ylab="spike alpha",
     xlim=c(0,1),ylim=yl,type="n",main="CT")
for (ii in 1:ncol(me.CT)) {
   lines(sim[,1],me.CT[,ii],col=ii)
}
palette(rainbow(ncol(me.CB)))
plot(0,0,xlab="Atwood number",ylab="spike alpha",
     xlim=c(0,1),ylim=yl,type="n",main="CB")
for (ii in 1:ncol(me.CB)) {
   lines(sim[,1],me.CB[,ii],col=ii)
}
palette(rainbow(ncol(me.CD)))
plot(0,0,xlab="Atwood number",ylab="spike alpha",
     xlim=c(0,1),ylim=yl,type="n",main="CD")
for (ii in 1:ncol(me.CD)) {
   lines(sim[,1],me.CD[,ii],col=ii)
}
graphics.off()

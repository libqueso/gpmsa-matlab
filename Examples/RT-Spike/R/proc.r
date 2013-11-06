sim <- as.matrix(read.table("../sim_spike.dat"))
obs <- as.matrix(read.table("../obs_spike.dat"))

obs <- obs[1:30,]
y.lim <- range(c(as.vector(sim[,-1]),obs[,2]))

postscript("spike.ps",height=4,width=4,paper="special",
           horizontal=F)
par(pty="s",mgp=c(2,1,0),oma=c(0,0,0,0),mar=c(3,3,1,1))
plot(0,0,xlab="Atwood number",ylab="spike alpha",xlim=c(0,1),
     ylim=y.lim,type="n")
for(i in 2:31) {
   lines(sim[,1],sim[,i],col="grey")
}
points(obs[,1],obs[,2],pch=19,col="cyan")
graphics.off()

PC <- sim[,2:31]
m.PC <- matrix(rep(apply(PC,1,mean),30),ncol=30)
PC <- PC - m.PC
sd.PC <- sd(as.vector(PC))
PC <- PC / sd.PC

svd.PC <- svd(PC)
pcd <- (svd.PC$d^2)/sum((svd.PC$d^2))
pcd <- round(100*pcd,2)
pct <- cumsum((svd.PC$d^2)/sum((svd.PC$d^2)))
# use 2 pc
K <- svd.PC$u[,1:2] %*% diag(svd.PC$d[1:2])
K.ns <- K*sd.PC
W <- svd.PC$v[,1:2]

pc1 <- matrix(0,nrow(K),3)
pc1[,1] <- m.PC[,1]+quantile(W[,1],.05)*K.ns[,1]
pc1[,2] <- m.PC[,1]
pc1[,3] <- m.PC[,1]+quantile(W[,1],.95)*K.ns[,1]

pc2 <- matrix(0,nrow(K),3)
pc2[,1] <- m.PC[,1]+quantile(W[,2],.05)*K.ns[,2]
pc2[,2] <- m.PC[,1]
pc2[,3] <- m.PC[,1]+quantile(W[,2],.95)*K.ns[,2]

postscript("spike_pc.ps",height=4,width=8,paper="special",horizontal=F)
par(mfrow=c(1,2),pty="s",mgp=c(2,1,0),oma=c(0,0,0,0),
    mar=c(3,3,1,1))
plot(sim[,1],pc1[,2],xlab="Atwood number",ylab="spike alpha",
     xlim=c(0,1),ylim=range(as.vector(pc1)),type="l",
     main=paste("PC 1 (", pcd[1], "%)",sep=""))
lines(sim[,1],pc1[,1],col="blue")
lines(sim[,1],pc1[,3],col="red")
plot(sim[,1],pc2[,2],xlab="Atwood number",ylab="spike alpha",
     xlim=c(0,1),ylim=range(as.vector(pc2)),type="l",
     main=paste("PC 2 (", pcd[2], "%)",sep=""))
lines(sim[,1],pc2[,1],col="blue")
lines(sim[,1],pc2[,3],col="red")
graphics.off()

kg <- seq(0,1,length=11)
kl <- matrix(0,nrow=nrow(obs),ncol=length(kg))
for (i in 1:ncol(kl)) {
   kl[,i] <- dnorm(obs[,1],mean=kg[i],sd=0.1)
}
m.obs <- approx(sim[,1],m.PC[,1],obs[,1])$y
vhat <- qr.solve(t(kl)%*%kl)%*%t(kl)%*%(obs[,2]-m.obs)
dhat <- kl%*%vhat

r = 0.01
R = 1+r
sig = 0.04
g4 = 1.01
beta = 10
g2 = -0.7
b2 = -0.4
g3 = 0.5
b3 = 0.3
g1 = 0
b1 = 0
b4 = 0

Y = matrix(NA, nrow=1, ncol=100)

for(i in 1:2){
  Y[i,1] = 0
  Y[i,2] = 0
  Y[i,3] = 0
  for(j in 4:100){
    U1 = (Y[i,j-1]-R*Y[i,j-2])*(g1*Y[i,j-3]+b1-R*Y[i,j-2])
    U2 = (Y[i,j-1]-R*Y[i,j-2])*(g2*Y[i,j-3]+b2-R*Y[i,j-2])
    U3 = (Y[i,j-1]-R*Y[i,j-2])*(g3*Y[i,j-3]+b3-R*Y[i,j-2])
    U4 = (Y[i,j-1]-R*Y[i,j-2])*(g4*Y[i,j-3]+b4-R*Y[i,j-2])
    
    l1 = beta*U1
    l2 =  beta*U2
    l3 = beta*U3
    l4 =  beta*U4
    
    N1 = 1+exp(l2-l1)+exp(l3-l1)+exp(l4-l1)
    N2 = exp(l1-l2)+1+exp(l3-l2)+exp(l4-l2) 
    N3 = exp(l1-l3)+exp(l2-l3)+1+exp(l4-l3)
    N4 = exp(l1-l4)+exp(l2-l4)+exp(l3-l4)+1
    
    n1 = 1/N1
    n2 = 1/N2
    n3 = 1/N3
    n4 = 1/N4
    
    Y[i,j] = 1/R*(n1*(g1*Y[i,j-1]+b1)+n2*(g2*Y[i,j-1]+b2)+n3*(g3*Y[i,j-1]+b3)+n4*(g4*Y[i,j-1]+b4)+rcauchy(1,scale=sig))
  }
}

Y_mc = Y

likelihood_log = function(theta, Y){
  likeli = numeric(length=100)
  likeli[1] = 0
  likeli[2] = 0
  likeli[3] = 0
  
  for(i in 4:100){
    p1 = exp(beta*(Y[i-1]-R*Y[i-2])*(g1*Y[i-3]+b1-R*Y[i-2]))
    p2 = exp(beta*(Y[i-1]-R*Y[i-2])*(theta[1]*Y[i-3]+theta[2]-R*Y[i-2]))
    p3 = exp(beta*(Y[i-1]-R*Y[i-2])*(theta[3]*Y[i-3]+theta[4]-R*Y[i-2]))
    p4 = exp(beta*(Y[i-1]-R*Y[i-2])*(g4*Y[i-3]+b4-R*Y[i-2]))
    p = p1+p2+p3+p4
    
    p11 = g1*Y[i-1]+b1
    p22 = theta[1]*Y[i-1]+theta[2]
    p33 = theta[3]*Y[i-1]+theta[4]
    p44 = g4*Y[i-1]+b4
    
    mu = (p1*p11+p2*p22+p3*p33+p4*p44)/(R*p)
    sd = sig/R
    
    likeli[i] = dnorm(Y[i], mu, sd, log=T)
  }
  likeli[is.na(likeli)]=0
  return(sum(likeli))
}

MH = function(x0, sd1, sd2, sd3, sd4, Niter){
  count = 0
  #posterior
  f = function(x) likelihood_log(theta=x, Y=Y_mc)+dunif(x[1], min=-2.5, max=0, log=T)+dunif(x[2], min=-2.5, max=0, log=T)+dunif(x[3], min=0, max=2.5, log=T)+dunif(x[4], min=0, max=1.5, log=T)
  
  #additive normal proposal
  rq = function(x) c(rnorm(1,x[1],sd1),rnorm(1,x[2],sd2), rnorm(1,x[3],sd3),rnorm(1,x[4],sd4))
  #exponential proposal
  #rq = function(x) c(-abs(x[1])*rexp(1,sd1), -abs(x[2])*rexp(1,sd2), abs(x[3])*rexp(1,sd3), abs(x[4])*rexp(1,sd4))
  
  #
  #q = function(x,y) dexp(y[1]/x[1], sd1)/x[1]*dexp(y[2]/x[2], sd2)/x[2]*dexp(y[3]/x[3], sd3)/x[3]*dexp(y[4]/x[4], sd4)/x[4]
  
  x = x0
  N = Niter
  xall = matrix(NA,nrow=4,ncol=N)
  xall[,1] = x
  for(i in 1:N-1){
    y = rq(x)
    alpha = exp(f(y)-f(x))
    #alpha = exp(f(y)-f(x))*(q(y,x)/q(x,y))
    
    if (runif(1)<alpha){
      count = count +1
      x = y
    }
    xall[,i+1] = x
  }
  burn = 2000
  #thin = 40
  xall = xall[,-c(1:burn)]
  #xall = xall[,seq(1,Niter-burn,by=thin)]
  return(xall)
}

X = MH(c(-0.93,-0.34,0.77,0.14),0.008,0.006,0.008,0.006,25000)
xall = X[,seq(1,25000-2000,by=50)]
write.table(xall, 'x11.csv', sep = ',')

pdf("trace1.pdf")


con1 = read.csv("C:/Users/dedaoyan/Desktop/xixi's project/con11.csv",head=T)
xall = read.csv("C:/Users/dedaoyan/Desktop/xixi's project/x11.csv",head=T)

par(mfrow=c(2,4))
plot(density(as.numeric(xall[1,])), ylab = "density", xlab = "g_2", col = "#56B4E9", xlim=c(-2.5,0),lwd = 2,  main=" ")
abline(v=-0.7, col="red", lwd = 2, lty=2)
abline(v=mean(as.numeric(xall[1,])), col="#56B4E9", lwd = 2, lty=2)
lines(x=seq(-2.5,0,length=10),rep(1/2.5,10), col="red", lwd = 2)
legend("topleft", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#56B4E9", "#56B4E9", "red", "red"), lty=c(1,2,1,2))

plot(density(as.numeric(xall[2,])), ylab = "density", xlab = "b_2", col = "#E69F00", lwd = 2, xlim=c(-1.5,0), main=" ")
abline(v=-0.4, col="red", lwd = 2, lty=2)
abline(v=mean(as.numeric(xall[2,])), col="#E69F00", lwd = 2, lty=2)
lines(x=seq(-1.5,0,length=10),rep(1/1.5,10), col="red", lwd = 2)
legend("topleft", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#E69F00", "#E69F00", "red", "red"), lty=c(1,2,1,2))

plot(density(con1$g2), ylab = "density", xlab = "g_2", col = "#56B4E9" ,lwd = 2, xlim=c(-4,1), main=" ")
abline(v=-0.7, col="red", lwd = 2, lty=2)
abline(v=mean(con1$g2), col="#56B4E9", lwd = 2, lty=2)
lines(x=seq(-2.5,0,length=10),rep(1/2.5,10), col="red", lwd = 2)
legend("topleft", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#56B4E9", "#56B4E9", "red", "red"), lty=c(1,2,1,2))

plot(density(con1$b2), ylab = "density", xlab = "b_2", col = "#E69F00", lwd = 2,xlim=c(-2,1),  main=" ")
abline(v=-0.4, col="red", lwd = 2, lty=2)
abline(v=mean(con1$b2), col="#E69F00", lwd = 2, lty=2)
lines(x=seq(-1.5,0,length=10),rep(1/1.5,10), col="red", lwd = 2)
legend("topleft", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#E69F00", "#E69F00", "red", "red"), lty=c(1,2,1,2))

plot(density(as.numeric(xall[3,])), ylab = "density", xlab = "g_3", col = "#009E73",xlim=c(0,2.5), lwd = 2,  main=" ")
abline(v=0.5, col="red", lwd = 3, lty=2)
abline(v=mean(as.numeric(xall[3,])), col="#009E73", lwd = 2, lty=2)
lines(x=seq(0,2.5,length=10),rep(1/2.5,10), col="red", lwd = 2)
legend("topright", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#009E73", "#009E73", "red", "red"), lty=c(1,2,1,2))

plot(density(as.numeric(xall[4,])), ylab = "density", xlab = "b_3", col = "#CC79A7",xlim=c(0,1.5), lwd = 2,  main=" ")
abline(v=0.3, col="red", lwd = 2, lty=2)
abline(v=mean(as.numeric(xall[4,])), col="#CC79A7", lwd = 2, lty=2)
lines(x=seq(0,1.5,length=10),rep(1/1.5,10), col="red", lwd = 2)
legend("topright", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#CC79A7", "#CC79A7", "red", "red"), lty=c(1,2,1,2))

plot(density(con1$g3), ylab = "density", xlab = "g_3", col = "#009E73", lwd = 2 ,xlim=c(-1,4), main=" ")
abline(v=0.5, col="red", lwd = 2, lty=2)
abline(v=mean(con1$g3), col="#009E73", lwd = 2, lty=2)
lines(x=seq(0,2.5,length=10),rep(1/2.5,10), col="red", lwd = 2)
legend("topright", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#009E73", "#009E73", "red", "red"), lty=c(1,2,1,2))

plot(density(con1$b3), ylab = "density", xlab = "b_3", col = "#CC79A7", lwd = 2, xlim=c(-0.1,3),main=" ")
abline(v=0.3, col="red", lwd = 2, lty=2)
abline(v=mean(con1$b3), col="#CC79A7", lwd = 2, lty=2)
lines(x=seq(0,1.5,length=10),rep(1/1.5,10), col="red", lwd = 2)
legend("topright", legend=c("posterior density", "posterior mean", "prior density", "true value"),
       col=c("#CC79A7", "#CC79A7", "red", "red"), lty=c(1,2,1,2))



mtext("            MH algorithm                                                                                      score based generative model", line=-3, cex=2, outer = TRUE)

par(mfrow = c(2,2))
plot(as.numeric(xall[1,]), type = "l", xlab = "iteration", ylab = "g_2", col = "#56B4E9", lwd = 2)
plot(as.numeric(xall[2,]), type = "l", xlab = "iteration", ylab = "b_2", col = "#E69F00", lwd = 2)
plot(as.numeric(xall[3,]), type = "l", xlab = "iteration", ylab = "g_3", col = "#009E73", lwd = 2)
plot(as.numeric(xall[4,]), type = "l", xlab = "iteration", ylab = "b_3", col = "#CC79A7", lwd = 2)
mtext("parameter set 1 with Cauchy noise", line=-3, cex=2, outer = TRUE)

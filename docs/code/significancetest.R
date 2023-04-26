rm(list=ls())

#===================
# HYPOTHESIS TESTS
#===================

#-----------------------------------
# power function
#-----------------------------------

mu_0 = 0
mus = seq(mu_0,mu_0+5,length.out = 100)
alpha = 0.1
n = 1
pow = pnorm(sqrt(n)*(mus - mu_0) - qnorm(1-alpha))

#pdf("Figure_power.pdf")
plot(mus, pow, type="l", ylim=c(0,1), xaxt="n", xlab="", yaxt="n", ylab=expression(pow(mu)))
axis(1,at=mu_0,labels=expression(mu[0]))
axis(2,at=c(0,alpha,1),labels=c(0,expression(alpha),1))
#dev.off()

#-----------------------------------
# interval hypothesis
#-----------------------------------

Delta = 0.5
n = 1:100
alpha = 0.1
c = sqrt( qchisq(alpha, df = 1, ncp = n*(Delta^2))/n )

#pdf("Figure_rejreg_interval.pdf")
plot(n,c, type="l", ylim=c(-1,1), 
     ylab=expression(bar(y)), yaxt="n")
axis(2, at=Delta, labels=expression(Delta))
axis(2, at=-Delta, labels=expression(-Delta))
axis(2, at=0)
lines(n,-c)
abline(h=-Delta, lty=2)
abline(h=Delta, lty=2)
#dev.off()

n = 25
c = sqrt( qchisq(alpha, df = 1, ncp = n*(Delta^2))/n )
mu = seq(-2,2,length.out = 101)
pow = 1-(pnorm(-c,mean=mu, sd=1/sqrt(n)) + pnorm(c,mean=mu, sd=1/sqrt(n), lower.tail = F))
#pdf("Figure_pow_interval.pdf")
plot(mu,pow, type="l", ylim=c(0,1), 
     xlab=expression(mu), ylab="Power", xaxt="n", yaxt="n")
axis(1, at=-Delta, labels=expression(-Delta))
axis(1, at=Delta, labels=expression(Delta))
axis(1, at=0, labels=0)
axis(2, at=alpha, labels=expression(alpha))
axis(2, at=c(0,1))
abline(h=alpha, lty=3)
#dev.off()

#-----------------------------------
# pulmonary data
#-----------------------------------

rm(list=ls())

library(ICSNP)
data(pulmonary)
Y = pulmonary
n = nrow(Y)
m = ncol(Y)
Y_bar = colMeans(Y)
S = var(Y)
t = n * t(Y_bar) %*% solve(S) %*% Y_bar
pf( ( (n-m)/(m*(n-1)) )* t, m,n-m , lower.tail = F)

id = c(2,3)
alpha = 0.05
Y = Y[,id]
n = nrow(Y)
m = ncol(Y)
Y_bar = colMeans(Y)
S = var(Y)

k <- sqrt((m/n)*((n-1)/(n-m)) * qf(1-alpha,m,n-m))
library(car)
#pdf("Figure_confreg.pdf")
plot(ellipse(center=Y_bar, shape=S, radius=k, draw=F), xlab=names(Y)[1], ylab=names(Y)[2], type="l")
ellipse(center=Y_bar, shape=S, radius=k, draw=T, col=1, cex=1)
abline(h=0, lty=3)
abline(v=0, lty=3)
#dev.off()

k_pred <- sqrt((m/n)*((n^2-1)/(n-m)) * qf(1-alpha,m,n-m))

#pdf("Figure_predreg.pdf")
plot(ellipse(center=Y_bar, shape=S, radius=k_pred, draw=F), xlab=names(Y)[1], ylab=names(Y)[2], type="l")
ellipse(center=Y_bar, shape=S, radius=k_pred, add=T, col=1, cex=1)
ellipse(center=Y_bar, shape=S, radius=k, add=T, col=1, cex=1, lty=2)
abline(h=0, lty=3)
abline(v=0, lty=3)
#dev.off()


Y = pulmonary
n = nrow(Y)
m = ncol(Y)
Y_bar = colMeans(Y)
S = var(Y)
c = qt(1-alpha/2, df=n-1)
d = sqrt( ( m*(n-1) / (n-m) ) * qf(1-alpha, m,n-m) )

U = Y_bar + diag(S)/sqrt(n) * c
L = Y_bar - diag(S)/sqrt(n) * c

U_s = Y_bar + diag(S)/sqrt(n) * d
L_s = Y_bar - diag(S)/sqrt(n) * d

round(cbind(Y_bar,L,U,L_s, U_s),2)


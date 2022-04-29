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

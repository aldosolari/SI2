rm(list=ls())

#===================
# SIGNIFICANCE TESTS
#===================

#-----------------------------------
# von Bortkiewicz's data
#-----------------------------------

#--- horse-kicks data 

y = c(rep(0,109),
      rep(1,65),
      rep(2,22),
      rep(3,3),
      4)

table(y)
n = length(y)

#-----------------------------------
# i) Test H_0: theta = theta_0
#-----------------------------------

theta_0 = 0.507
alpha = 0.05

#--- log-likelihood

x <- sum(y)
lambda_0 <- n * theta_0 
lambda_hat <- x

loglik <- function(lambda) {
  ll <- x * log(lambda) - lambda
  return(ll)
}

a <- 70
b <- 150
#pdf("Figure_loglikelihood.pdf")
curve(loglik,
      from=a,to=b, 
      xlab=expression(lambda), 
      ylab=expression(loglikelihood(lambda)),
      xaxt="n", yaxt="n",
      lwd=2)
axis(1,at=c(a,lambda_0,lambda_hat,b), 
     labels=c(a,expression(lambda[0]), expression(hat(lambda)),b))
axis(2,at=c(loglik(lambda_0),loglik(lambda_hat)), labels=c("",""))
segments(x0=lambda_hat,x1=lambda_hat,y0=0,y1=loglik(lambda_hat), lty=3)
segments(x0=lambda_0,x1=lambda_0,y0=0,y1=loglik(lambda_0), lty=3)
segments(x0=0,x1=lambda_0,y0=loglik(lambda_0),y1=loglik(lambda_0), lty=3)
segments(x0=0,x1=lambda_hat,y0=loglik(lambda_hat),y1=loglik(lambda_hat), lty=3)
abline(b=(lambda_hat-lambda_0)/(lambda_0), a=loglik(lambda_0)-((lambda_hat-lambda_0)/(lambda_0))*lambda_0, lty=3)
#dev.off()

#--- likelihood-based tests and CI

# wald
t_E = (lambda_hat - lambda_0) / sqrt(lambda_hat)
p_wald = 2*pnorm(abs(t_E), lower.tail = F)

z_alpha = -qnorm(alpha/2)
ci_wald = (lambda_hat + c(-1,1)*z_alpha*sqrt(lambda_hat) ) 
ci_wald / n

# score
t_U = (lambda_hat - lambda_0) / sqrt(lambda_0)
p_score = 2*pnorm(abs(t_U), lower.tail = F)

ci_score = (lambda_hat + z_alpha^2 / 2 + c(-1, 1) * z_alpha * sqrt(lambda_hat + z_alpha^2 / 4) )
ci_score / n

# likelihood ratio
t_L = 2 * (loglik(lambda_hat) - loglik(lambda_0))
p_lrt = pchisq(t_L, lower.tail = F, df=1)

c_alpha = qchisq(alpha, df = 1, lower.tail = F)
ci_lrt <- vector()
ci_lrt[1]<-uniroot(function(x) loglik(lambda_hat) - loglik(x) - c_alpha/2, lower = 0, upper = lambda_hat)$root 
ci_lrt[2]<-uniroot(function(x) loglik(lambda_hat) - loglik(x) - c_alpha/2, lower = lambda_hat, upper = n)$root 
ci_lrt / n

exp(confint(glm(y ~ 1, family=poisson)))

a <- 90
b <- 150
#pdf("Figure_lr_interval.pdf")
curve(loglik,
      from=a,to=b, 
      xlab=expression(lambda), 
      ylab=expression(loglikelihood(lambda)),
      xaxt="n", yaxt="n",
      lwd=2)
axis(1,at=c(ci_lrt[1],lambda_0,lambda_hat,ci_lrt[2]), 
     labels=c("",expression(lambda[0]), expression(hat(lambda)), ""))
axis(2,at=loglik(lambda_hat) - c_alpha/2, label=expression(l(hat(lambda))-c[alpha]/2))
abline(h = loglik(lambda_hat) - c_alpha/2, lty=2)
segments(x0=ci_lrt[1],x1=ci_lrt[1],y0=0,y1=loglik(ci_lrt[1]), lty=3)
segments(x0=ci_lrt[2],x1=ci_lrt[2],y0=0,y1=loglik(ci_lrt[1]), lty=3)
#dev.off()

#--- poisson.test

t = sum(y)
poisson.test(x=t, T=n, r=theta_0, conf.level = 1-alpha)

#--- max(Y_i) test

library(expint)
1 - ( gammainc(max(y), x=theta_0)^n ) / ( gamma(max(y)) ^ n )

#--- two-sided test

p_plus_obs = 1-ppois(t-1, lambda = lambda_0)
p_minus_obs = ppois(t, lambda = lambda_0)

other_tail <- ppois( qpois(p_plus_obs, lambda_0) - 1, lambda = lambda_0)
q_obs = min(p_plus_obs,p_plus_obs) + other_tail
q_obs

#--- confidence interval for a Poisson mean

p_theta_garwood = function(theta,t){
  p_plus <- ppois(t - 1, theta, lower.tail = FALSE)
  p_minus <- ppois(t, theta, lower.tail = TRUE)
  p = min(2*min(p_plus,p_minus),1)
}

p_theta_blaker = function(theta,t){
   if (theta <= t) {
     p_plus <- ppois(t - 1, theta, lower.tail = FALSE)
     other_tail <- ppois( qpois(p_plus, theta) - 1, theta)
     p <- min(p_plus + other_tail, 1)
   }
   else {
     p_minus <- ppois(t, theta, lower.tail = TRUE)
     other_tail <- ppois( qpois(1-p_minus, theta), theta, lower.tail = FALSE) 
     p <- min(p_minus + other_tail, 1)
   }
return(p)
}

#--- garwood

qchisq(alpha/2, df=2*t) / (2*n)
qchisq(1-alpha/2, df=2*(t+1)) / (2*n)

uniroot(function(x) p_theta_garwood(theta=x,t=t) - alpha, interval = c(0,t))$root / n
uniroot(function(x) p_theta_garwood(theta=x,t=t) - alpha, interval = c(t,n))$root / n

#--- blaker

uniroot(function(x) p_theta_blaker(theta=x,t=t) - alpha, interval = c(0,t))$root / n
uniroot(function(x) p_theta_blaker(theta=x,t=t) - alpha, interval = c(t,n))$root / n

theta_seq = seq(0.4,0.8,length.out = 10^4)
p_garwood = sapply(theta_seq, function(theta) p_theta_garwood(n*theta,t = t))
p_blaker = sapply(theta_seq, function(theta) p_theta_blaker(n*theta,t = t))

#pdf("Figure_pvalue_function.pdf")
plot(theta_seq, p_garwood, 
     type="s", ylim=c(0,1), col="gray",
     xlab=expression(theta), ylab=expression(p[theta]))
lines(theta_seq, p_blaker)
abline(h=alpha, lty=2)
legend("topright", c("Garwood","Blaker"), col=c("gray",1), lty=1)
#dev.off()

#--- coverage

theta_seq = seq(0,10,length.out = 10^4)
t_seq = 0:20

lo_garwood = qchisq(alpha/2, df=2*t_seq)/2
up_garwood = qchisq(1-alpha/2, df=2*(t_seq+1)) / 2

library(BlakerCI)
ci_blaker = sapply(t_seq, function(x) poisson.blaker.limits(x, level = 1-alpha))

coverage_garwood = sapply(theta_seq, function(theta) 
  sum(dpois(t_seq[theta >= lo_garwood & theta <= up_garwood], lambda=theta)))

coverage_blaker = sapply(theta_seq, function(theta) 
  sum(dpois(t_seq[theta >= ci_blaker[1,] & theta <= ci_blaker[2,]], lambda=theta)))

#pdf("Figure_intervals.pdf")
plot(t_seq, lo_garwood, type="s", xlab = expression(t[obs]), ylab="Confidence interval", col=0, ylim=c(min(lo_garwood),max(up_garwood)) )
for (i in 1:length(t_seq)) segments(x0=t_seq[i],x1=t_seq[i],y0=lo_garwood[i], y1=up_garwood[i], col="gray")
for (i in 1:length(t_seq)) segments(x0=t_seq[i]+.2,x1=t_seq[i]+.2,y0=ci_blaker[1,i], y1=ci_blaker[2,i])
legend("topleft", c("Garwood","Blaker"), col=c("gray",1), lty=1)
#dev.off()

#pdf("Figure_coverage.pdf")
plot(theta_seq, coverage_garwood, 
     type="s", col="gray",
     ylim = c(1-alpha,1),
     xlab=expression(theta), ylab="Coverage probability")
lines(theta_seq, coverage_blaker, type="s")
abline(h=1-alpha, lty=2)
legend("topright", c("Garwood","Blaker"), col=c("gray",1), lty=1)
#dev.off()

#-----------------------------------
# ii) Test H_0: Y_i ~ Poisson(theta)
#-----------------------------------

s = sum(y)

set.seed(123)
B = 10^4 - 1
t_dispersion <- vector()
t_zeros <- vector()
for (b in 1:B){
  y_new = rmultinom(n=1, size=s, prob=rep(1/n,n))
  t_dispersion[b] = (n-1)*var(y_new) / mean(y_new)
  t_zeros[b] = sum(y_new==0)
}
t_dispersion[B+1] <- (n-1)*var(y) / mean(y)
t_zeros[B+1] <- sum(y==0)

#pdf("Figure_dispersion.pdf")
plot(table(t_dispersion), xaxt="n", xlab="Dispersion index", ylab="Frequency")
axis(1,at=t_dispersion[B+1], labels=expression(t[obs]))
axis(1,at=range(t_dispersion), labels=round(range(t_dispersion),1))
#dev.off()

mean(t_dispersion >= t_dispersion[B+1])
pchisq(t_dispersion[B+1], df=n-1, lower.tail = F)

#pdf("Figure_zeros.pdf")
plot(table(t_zeros), xaxt="n", xlab="Number of zeros", ylab="Frequency")
axis(1,at=t_zeros[B+1], labels=expression(t[obs]))
axis(1,at=range(t_zeros), labels=round(range(t_zeros),1))
#dev.off()

mean(t_zeros >= t_zeros[B+1])

#-----------------------------------
# Darwin’s data
#-----------------------------------

library(scifigure)
#pdf("Figure_sci_exp.pdf")
sci_figure(init_experiments(1), showlegend = F)
#dev.off()

library(SMPracticals)
library(ggpubr)

#--- two-sample t test

darwin_pair = data.frame(Pot = darwin[darwin$type=="Cross",1],
                         Cross = darwin[darwin$type=="Cross",4],
                         Self = darwin[darwin$type=="Self",4])
darwin_pair

#--- two-sample t test

#pdf("Figure_darwin_boxplot.pdf")
ggboxplot(darwin, x = "type", y = "height",color = "type")
#dev.off()

t.test(height  ~ type, data=darwin, var.equal=TRUE)

#--- paired t test

#pdf("Figure_darwin_diff.pdf")
ggline(darwin, x = "type", y = "height", group="pair")
#dev.off()

z <- apply(darwin_pair[,3:2],1,diff)
t.test(z)

#--- test of symmetry

n = length(z)
signs <- e1071::bincombinations(n)
signs[which(signs == 1)] <- -1
signs[which(signs == 0)] <- 1

ts = sapply(1:nrow(signs), function(i) z %*% signs[i,])

#pdf("Figure_symmetry.pdf")
plot(table(ts), xaxt="n", ylab="Frequency", xlab="")
axis(1,at=ts[1], labels=expression(t[obs]))
#dev.off()

p_plus = mean(ts >= ts[1])
p_minus = mean(ts <= ts[1])
2*min(p_plus, p_minus)

#--- Wilcoxon signed-rank test

z_abs = abs(z)
sr = sum(rank(z_abs) * sign(z))

ws = sapply(1:nrow(signs), function(i) sum( rank(z_abs) * sign(z) * signs[i,]))
#pdf("Figure_signedrank.pdf")
plot(table(ws), xaxt="n", ylab="Frequency", xlab="")
axis(1,at=ws[1], labels=expression(t[obs]))
#dev.off()

p_plus = mean(ws >= ws[1])
p_minus = mean(ws <= ws[1])
2*min(p_plus, p_minus)
wilcox.test(darwin_pair$Cross, darwin_pair$Self, paired=TRUE) 

#--- Sign test

#pdf("Figure_sign.pdf")
plot(0:n, dbinom(0:n,size=n,p=0.5), type="h", ylab="Probability", xlab="")
axis(1,at=sum(z>0), labels=expression(t[obs]))
#dev.off()

p_minus = pbinom(sum(z>0), size=n, p=0.5)
p_plus = pbinom(sum(z>0)-1, size=n, p=0.5, lower.tail = F)
2*min(p_plus, p_minus)

binom.test(x=sum(z>0), n=length(z), p=0.5)



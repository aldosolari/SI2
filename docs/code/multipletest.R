rm(list=ls())

#===================
# MULTIPLE TESTING
#===================

#-----------------------------------
# p.adjust example
#-----------------------------------

set.seed(123)
x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
p <- 2*pnorm(sort(-abs(x)))

round( p.adjust(p,"bonferroni"), 3)


#-----------------------------------
# Simulation
#-----------------------------------

rm(list=ls())

sim_errors <- function(m, m0, effect, method="none") {
  alpha = 0.05
  stats <- rnorm(m)
  stats[(m0+1):m] <- stats[(m0+1):m] + effect
  pvals = pnorm(stats, lower.tail = FALSE)
  setR <- which(p.adjust(pvals,method) <= alpha)
  # number of rejections
  R = length(setR)
  # number of type I errors
  V <- sum( setR <= m0 )
  # family-wise error
  FWER <- (V > 0)
  # false discovery proportion
  FDP <- V/max(1,R)
  return(c(R, V, FWER, FDP))
}

set.seed(123)
B = 10^5
res <- replicate(B, sim_errors(20, 15, 2, method="none"))
row.names(res) <- c("R","V","V>0", "V/R")
res[,1:10]

xtable::xtable(res[,1:10])

# FWER
mean(res[3,])
# FDR 
mean(res[4,])



#-----------------------------------
# fwer control
#-----------------------------------

rm(list=ls())

m0 = 100
FWER = 1-(1-0.05)^(1:m0)
#pdf("Figure_fwer.pdf")
plot(1:m0,FWER, type="l", xlab=expression(m[0]), ylab="FWER", lwd=2)
#dev.off()



######################################
## Multiple comparison with control ##
######################################

## (Dunnnet)

library(mvtnorm) ## for multivariate t-distribution

qDunnett <- function(q, g, df, tail){
  ## q:    probability
  ## g:    number of groups
  ## tail: see ?qmvt  
  
  corr       <- matrix(0.5, nrow = g - 1, ncol = g - 1) ## implicitly assumes *balanced* design!
  diag(corr) <- 1
  
  if(length(g) > 1)
    stop("only scalar g allowed")
  
  nq  <- length(q)
  out <- numeric(nq)
  
  
  for(j in seq_along(q)){
    out[j] <- qmvt(q, tail = tail, df = df, corr = corr)$quantile
    ##uniroot(function(x, g) eval_t(x, g) - q[j], g = g, 
    ##interval = c(0, 20))$root
  }
  out
}

g.sizes <- c(3:11, 16)#, 21, 31, 41)

level <- 0.95

critval.two <- critval.one <- numeric(length(g.sizes))

for(j in seq_along(g.sizes)){
  critval.two[j] <- qDunnett(level, g.sizes[j], df = 10, tail = "both.tails") 
  critval.one[j] <- qDunnett(level, g.sizes[j], df = 10, tail = "lower.tail") 
}

round(critval.two, 2) ## compare with Oehlert, page 637 => looks good, only minor differences
round(critval.one, 2) ## compare with Oehlert, page 635 => looks good, only minor differences


## Other package with quantile functions ####
library(nCDunnett)
quantiles <- sapply(g.sizes, function(x) qNCDun(0.95, nu = 10, rho = rep(0.5, x - 1), delta = 0, n = 32))
round(quantiles, 2)

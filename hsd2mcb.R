## See Hsu (1996), page. 108 for an example or Theorem 4.2.1 on page 104 for the theory.

## Method based on pairwise comparison a la Tukey HSD ####
hsd2mcb <- function(x){ 
  nx  <- length(x)
  out <- setNames(vector("list", nx), names(x))
  
  ## needs some more careful error checking for strange level names containing "-"
  tmp      <- attr(x, "orig.call")
  tmp[[1]] <- as.name("model.frame")
  mf       <- eval(tmp)
  
  if(any(attr(attr(mf, "terms"), "order") > 1))
    stop("Currently only one-way ANOVA models supported")

  for(j in names(x)){
    mat   <- x[[j]]
    nms   <- levels(mf[,j])
    ci    <- matrix(0, nrow = length(nms), ncol = 2)
    colnames(ci) <- c("lwr", "upr")
    rownames(ci) <- paste(nms, "-max(other)", sep = "")
    d     <- rownames(mat)
    for(i in seq_along(nms)){
      left  <- grep(paste("^", nms[i], "-", sep = ""), d)
      right <- grep(paste("-", nms[i], "$", sep = ""), d)
      ci[i,] <- c(min(mat[left ,"lwr"], -mat[right, "upr"]),
                  min(mat[left ,"upr"], -mat[right, "lwr"]))
    }
    out[[j]] <- ci
  }
  out
}

## Data from Hsu, analysis on page 107
data <- data.frame(y = c(45, 59, 48, 46, 38, 47,
                         21, 12, 14, 17, 13, 17,
                         37, 32, 15, 25, 39, 41,
                         16, 11, 20, 21, 14, 7),
                   color = rep(c("yellow", "white", "red", "blue"), each = 6))
fit <- aov(y ~ color, data = data)
hsd <- TukeyHSD(fit, conf.level = 0.99)
hsd2mcb(hsd)

## R data-set with *two* variables, see ?TukeyHSD
summary(fm1 <- aov(breaks ~ wool + tension, data = warpbreaks))
fm1.hsd <- TukeyHSD(fm1, c("wool", "tension"))
hsd2mcb(fm1.hsd)

## including interaction
fm2 <- aov(breaks ~ wool * tension, data = warpbreaks)
fm2.hsd <- TukeyHSD(fm2)
hsd2mcb(fm2.hsd)
## error because of interaction term (as intended)

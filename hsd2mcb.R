## See Hsu (1996), page. 108 for an example or Theorem 4.2.1 on page 104 for the theory.

## Method based on pairwise comparison a la Tukey HSD ####

# using TukeyHSD object
hsd2mcb.tuk <- function(x){ 
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

# using glht object
hsd2mcb <- function(x, conf.level = .95){ 
  
  nx       <- length(x$focus)
  factor_names <- x$focus
  out      <- setNames(vector("list", nx), factor_names)
  
  ## needs some more careful error checking for strange level names containing "-"
  
  tmp <- class(x$model)
  
  if(!any(tmp %in% c("aov","lm","lmerMod","lme4","glm")))
    stop("Currently only aov, lm, glm, lmerMod are supported")
  
  mf <- if(any(tmp %in% c("lmerMod", "lme4"))) x$model@frame else x$model$model

  if(any(attr(attr(mf, "terms"), "order") > 1))
    stop("Currently no interactions are supported")
  
  # confindence intervals (ci)
  mat_temp <- confint(x, level = conf.level)$confint
  mat_all <- list(mat_temp)
  
  # if testing more than one factor, split ci matrix by factor
  if(length(factor_names)>1){
    
    # get factor names
    fnames <- gsub(pattern = "^(.+)(: )(.+)$", 
                   replacement = "\\1", 
                   x = rownames(mat_temp), perl = TRUE)
    # leave only level names
    cnames <- gsub(pattern = "^(.+: )", 
                   replacement = "", 
                   x = rownames(mat_temp), perl = TRUE)
    rownames(mat_temp) <- cnames
    
    mat_all <- list()
    
    # assign a matrix per factor
    for(n in unique(fnames)){
      mat_all[[n]] <- mat_temp[fnames %in% n, ,drop=FALSE]
    }
    
  }
  
  names(mat_all) <- factor_names
  
  for(j in factor_names){
    
    mat <- mat_all[[j]]
    nms <- levels(mf[,j])
    ci <- matrix(0, nrow = length(nms), ncol = 2)
    colnames(ci) <- c("lwr", "upr")
    rownames(ci) <- paste(nms, "-max(other)", sep = "")
    d <- rownames(mat)
    
    for(i in seq_along(nms)){
      
      left  <- grep(paste("^", nms[i], " - ", sep = ""), d)
      right <- grep(paste(" - ", nms[i], "$", sep = ""), d)
      
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
fit.hsd <- TukeyHSD(fit, conf.level = 0.99)
library(multcomp)
fit.glht <- glht(fit, linfct = mcp(color = "Tukey"))
hsd2mcb.tuk(fit.hsd)
hsd2mcb(fit.glht, conf.level = .99)


## R data-set with *two* variables, see ?TukeyHSD
summary(fm1 <- aov(breaks ~ wool + tension, data = warpbreaks))
fm1.hsd <- TukeyHSD(fm1, c("wool", "tension"))
fm1.glht <- glht(fm1, linfct = mcp(wool = "Tukey", tension = "Tukey"))
fm1.glhtw <- glht(fm1, linfct = mcp(wool = "Tukey"))
fm1.glhtt <- glht(fm1, linfct = mcp(tension = "Tukey"))
hsd2mcb.tuk(fm1.hsd)
hsd2mcb(fm1.glht) # not the same because of simultaneous confidence intervals
hsd2mcb(fm1.glhtw)
hsd2mcb(fm1.glhtt)

## including interaction
#+ error = TRUE
fm2 <- aov(breaks ~ wool * tension, data = warpbreaks)
fm2.hsd <- TukeyHSD(fm2)
fm2.glht <- glht(fm2, linfct = mcp(wool = "Tukey", tension = "Tukey"))
hsd2mcb.tuk(fm2.hsd)
hsd2mcb(fm2.glht)
## error because of interaction term (as intended)

warpbreaks$wooltension <- factor(paste(warpbreaks$wool, warpbreaks$tension))
fm3 <- aov(breaks ~ wooltension, data = warpbreaks)
fm3.hsd <- TukeyHSD(fm3)
fm3.glht <- glht(fm3, linfct = mcp(wooltension = "Tukey"))
hsd2mcb.tuk(fm3.hsd)
hsd2mcb(fm3.glht)

## glm example from help page
utils::data(anorexia, package = "MASS")

anorex.1 <- glm(Postwt ~ Prewt + Treat + offset(Prewt),
                family = gaussian, data = anorexia)
summary(anorex.1) # FT has largest mean

glm.glht <- glht(anorex.1, linfct = mcp(Treat = "Tukey"))
mat <- confint(glm.glht)$confint
true <- rbind(c(min(-mat["Cont - CBT","upr"], -mat["FT - CBT","upr"]),
                min(-mat["Cont - CBT","lwr"], -mat["FT - CBT","lwr"])),
              c(min(mat["Cont - CBT","lwr"], -mat["FT - Cont","upr"]),
                min(mat["Cont - CBT","upr"], -mat["FT - Cont","lwr"])),
              c(min(mat["FT - CBT","lwr"], mat["FT - Cont", "lwr"]),
                min(mat["FT - CBT","upr"], mat["FT - Cont", "upr"])))

dimnames(true) <- list(c("CBT - max(other)", 
                         "Cont - max(other)", 
                         "FT - max(other)"),
                       c("lwr","upr"))
true
hsd2mcb(glm.glht)

## lmer
library(lme4)

data(Machines,package="nlme")

mach <- lmer(score ~ Machine + (Machine-1|Worker), data=Machines)
summary(mach) # would expect Machine C to be the max

lmer.glht <- glht(mach, linfct = mcp(Machine = "Tukey"))
hsd2mcb(lmer.glht)

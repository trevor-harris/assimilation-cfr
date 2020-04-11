#### QI

# Unvariate Tukey Depth
depth = function(g, fmat) {
  
  # Computes the depth values of a function with respect to a set of functions (fmat)
  fn = ncol(fmat)
  depth = rep(0, length(g))
  
  for (row in 1:nrow(fmat)) {
    diff = abs(sum(sign(g[row] - fmat[row,])))
    depth[row] = 1 - (diff / fn)
  }
  
  return(depth)
}

# Integrated Tukey Depth
xdepth = function(gmat, fmat) {
  apply(gmat, 2, function(x) mean(depth(x, fmat)))
}

# Quality Index with Tukey Depth
# f and g are two sets of functions, represented as matricies where each column is a separate function 
quality = function(f, g) {
  fxd = xdepth(f, f)
  gxd = xdepth(g, f)
  
  q = mean(sapply(gxd, function(y) mean(fxd <= y)))
  pval = 1 - pnorm((0.5 - q) / sqrt(((1/ncol(f)) + (1/ncol(g)))/12))
  
  c(q, pval)
}


#### FAD
adk.test <-
  function (...) 
  {
    #############################################################################
    # This function "adk.test" tests whether k samples (k>1) come from a common
    # continuous distribution, using the nonparametric (rank) test described in
    # Scholz F.W. and Stephens M.A. (1987), K-sample Anderson-Darling Tests,
    # Journal of the American Statistical Association, Vol 82, No. 399, 
    # pp. 918-924. 
    # This test is consistent against all alternatives. 
    # There is an adjustment for ties, and for a moderate number of tied 
    # observations the P-value calculation should be reasonable.
    #
    # When there are NA's among the sample values they are removed,
    # with a warning message indicating the number of NA's.
    # It is up to the user to judge whether such removals make sense.
    # 
    # The inputs can either be a sequence or a list of k (>1) sample vectors.
    #
    # An example:
    # x1 <- c(1, 3, 2, 5, 7), x2 <- c(2, 8, 1, 6, 9, 4), 
    # and x3 <- c(12,  5,  7,  9, 11)
    # adk.test(x1,x2,x3) # or adk.test(list(x1,x2,x3)) produces the output below.
    #############################################################################
    # Anderson-Darling k-sample test.
    #
    # Number of samples:  3
    # Sample sizes: 5 6 5
    # Total number of values: 16
    # Number of unique values: 11
    #
    # Mean of Anderson-Darling Criterion: 2
    # Standard deviation of Anderson-Darling Criterion: 0.92837
    #
    # T.AD = (Anderson-Darling Criterion - mean)/sigma
    #
    # Null Hypothesis: All samples come from a common population.
    #
    #                      T.AD P-value extrapolation
    # not adj. for ties 1.41756 0.08956             0
    # adj. for ties     1.62856 0.07084             0
    #############################################################################
    # In order to get the output list, call out.adk <- adk.test(x1,x2,x3)
    # then out.adk is of class adk and has components 
    # > names(out.adk)
    # [1] "k"       "ns"      "n"       "n.ties"  "sig"    "adk"     "warning"
    #
    # where
    # k = number of samples being compared
    # ns = vector of the k sample sizes
    # n = total sample size = n_1+...+n_k
    # n.ties = number of ties in the combined set of all n observations
    # sig = standard deviation of the AD statistic
    # adk = 2 x 3 matrix containing t.obs, P-value, extrapolation, adjusting for
    #       ties and not adjusting for ties. extrapolation is TRUE when 
    #       the P-value was extrapolated. Here t.obs is the observed value of 
    #       T.AD and P-value = P(T.AD > t.obs).
    # warning = logical indicator, warning = TRUE indicates that at least  
    # one of the sample sizes is < 5.
    #
    # Fritz Scholz, May 2011
    #
    #################################################################################
    na.remove <- function(x){
      #
      # This function removes NAs from a list and counts the total 
      # number of NAs in na.total.
      # Returned is a list with the cleaned list x.new and with 
      # the count na.total of NAs.
      #
      na.status <- sapply(x,is.na)
      k <- length(x)
      x.new <- list()
      na.total <- 0
      for( i in 1:k ){
        x.new[[i]] <- x[[i]][!na.status[[i]]]
        na.total <- na.total + sum(na.status[[i]])
      }
      list(x.new=x.new,na.total=na.total)
    }
    if (nargs() == 1 & is.list(list(...)[[1]])) {
      samples <- list(...)[[1]]
    }
    else {
      samples <- list(...)
    }
    out <- na.remove(samples)
    na.t <- out$na.total
    if( na.t > 1) cat(paste("\n",na.t," NAs were removed!\n\n"))
    if( na.t == 1) cat(paste("\n",na.t," NA was removed!\n\n"))
    samples <- out$x.new
    k <- length(samples)
    if (k < 2) 
      stop("Must have at least two samples.")
    ns <- sapply(samples, length)
    if (any(ns == 0)) 
      stop("One or more samples have no observations.")
    x <- NULL
    for (i in 1:k) x <- c(x, samples[[i]])
    n <- length(x)
    Z.star <- sort(unique(x))
    L <- length(Z.star)
    AkN2 <- 0
    AakN2 <- 0
    l.vec <- NULL
    for (j in 1:L) {
      fij <- NULL
      for (i in 1:k) {
        fij <- c(fij, as.double(sum(samples[[i]] == Z.star[j])))
      }
      l.vec <- c(l.vec, sum(fij))
    }
    for (i in 1:k) {
      Mij <- as.double(0)
      Maij <- as.double(0)
      inner.sum <- as.double(0)
      inner.sum.a <- as.double(0)
      for (j in 1:L) {
        fij <- sum(samples[[i]] == Z.star[j])
        Mij <- Mij + fij
        Maij <- Mij - fij/2
        Bj <- as.double(sum(l.vec[1:j]))
        Baj <- Bj - l.vec[j]/2
        if (j < L) {
          inner.sum <- inner.sum + l.vec[j] * (n * Mij - 
                                                 ns[i] * Bj)^2/(Bj * (n - Bj))
        }
        inner.sum.a <- inner.sum.a + l.vec[j] * (n * Maij - 
                                                   ns[i] * Baj)^2/(Baj * (n - Baj) - n * l.vec[j]/4)
      }
      AkN2 <- AkN2 + inner.sum/ns[i]
      AakN2 <- AakN2 + inner.sum.a/ns[i]
    }
    AkN2 <- AkN2/n
    AakN2 <- (n - 1) * AakN2/n^2
    coef.d <- 0
    coef.c <- 0
    coef.b <- 0
    coef.a <- 0
    H <- sum(1/ns)
    h <- sum(1/(1:(n - 1)))
    g <- 0
    for (i in 1:(n - 2)) {
      g <- g + (1/(n - i)) * sum(1/((i + 1):(n - 1)))
    }
    coef.a <- (4 * g - 6) * (k - 1) + (10 - 6 * g) * H
    coef.b <- (2 * g - 4) * k^2 + 8 * h * k + (2 * g - 14 * h - 
                                                 4) * H - 8 * h + 4 * g - 6
    coef.c <- (6 * h + 2 * g - 2) * k^2 + (4 * h - 4 * g + 6) * 
      k + (2 * h - 6) * H + 4 * h
    coef.d <- (2 * h + 6) * k^2 - 4 * h * k
    sig2 <- (coef.a * n^3 + coef.b * n^2 + coef.c * n + coef.d)/((n - 
                                                                    1) * (n - 2) * (n - 3))
    sig <- sqrt(sig2)
    TkN <- (AkN2 - (k - 1))/sig
    TakN <- (AakN2 - (k - 1))/sig
    pvalTkN <- adk.pval(TkN, k - 1)
    pvalTakN <- adk.pval(TakN, k - 1)
    warning <- min(ns) < 5
    ad.mat <- matrix(c(signif(TkN, 7), round(pvalTkN[[1]][1], 
                                             7), pvalTkN[[2]][1], signif(TakN, 7), round(pvalTakN[[1]][1], 
                                                                                         7), pvalTakN[[2]][1]), byrow = TRUE, ncol = 3)
    dimnames(ad.mat) <- list(c("not adj. for ties", "adj. for ties"), 
                             c("t.obs", "P-value", "extrapolation"))
    object <- list(k = k, ns = ns, n = n, n.ties = n - L, sig = round(sig, 
                                                                      5), adk = round(ad.mat, 5), warning = warning)
    class(object) <- "adk"
    object
  }
adk.pval <-
  function (tx,m) 
  {
    # This function "adk.pval" evaluates the p-value of the observed value 
    # tx of the T_m statistic in "K-Sample Anderson-Darling Tests" by F.W. Scholz 
    # and M.A. Stephens (1987), Journal of the American Statistical Association, 
    # Vol 82, No. 399, pp 918-924. Thus this p-value is P(T_m >= tx).
    #
    # Input: tx = observed value of T_m, tx > 0.
    #         m = the index of T_m, m >= 1.
    #
    # Output: a list with components
    #         p0 = p-value of tx, i.e., p0 = P(T_m >= tx)
    #         extrap = a logical indicator
    #                    extrap = TRUE means that linear extrapolation took place
    #                    extrap = FALSE means that quadratic interpolation was used.
    #
    # Computational Details:
    #
    # This function first interpolates the upper T_m quantiles as given in Table 1
    # of the above reference to the given value of m by fitting a quadratic 
    # in 1/sqrt(m) to the quantiles as tabulated for the upper quantile
    # levels .25, .10, .05, .025, .01.
    #
    # Next a quadratic in the interpolated quantiles (for m) is fitted to the  
    # log-odds of the upper probability levels defining these quantiles
    # and the fitted log-odds value at tx is converted back to the calculated upper 
    # probability value, i.e., the p-value. p-values outside the tabulated range
    # [.01,.25] are obtained by linear extrapolation of the fitted quadratic.
    #  
    # Fritz Scholz, May 2011
    #=========================================================================================
    table1.adk <- cbind(c(1, 2, 3, 4, 6, 8, 10, Inf), c(0.326, 
                                                        0.449, 0.498, 0.525, 0.557, 0.576, 0.59, 0.674), c(1.225, 
                                                                                                           1.309, 1.324, 1.329, 1.332, 1.33, 1.329, 1.282), c(1.96, 
                                                                                                                                                              1.945, 1.915, 1.894, 1.859, 1.839, 1.823, 1.645), c(2.719, 
                                                                                                                                                                                                                  2.576, 2.493, 2.438, 2.365, 2.318, 2.284, 1.96), c(3.752, 
                                                                                                                                                                                                                                                                     3.414, 3.246, 3.139, 3.005, 2.92, 2.862, 2.326))
    extrap <- FALSE
    mt <- table1.adk[, 1]
    sqm1 <- 1/sqrt(mt)
    sqm2 <- sqm1^2
    tm <- NULL
    p <- c(0.25, 0.1, 0.05, 0.025, 0.01)
    lp <- log(p/(1 - p))
    for (i in 1:5) {
      out <- lsfit(cbind(sqm1, sqm2), table1.adk[, i + 1])
      x <- 1/sqrt(m)
      coef <- out$coef
      y <- coef[1] + coef[2] * x + coef[3] * x^2
      tm <- c(tm, y)
    }
    out <- lsfit(cbind(tm, tm^2), lp)
    coef <- out$coef
    if (tx <= max(tm) & tx >= min(tm)) {
      lp0 <- coef[1] + coef[2] * tx + coef[3] * tx^2
    }
    if (tx > max(tm)) {
      extrap <- TRUE
      lp0 <- min(lp) + (tx - max(tm)) * (coef[2] + 2 * coef[3] * 
                                           max(tm))
    }
    if (tx < min(tm)) {
      extrap <- TRUE
      lp0 <- max(lp) + (tx - min(tm)) * (coef[2] + 2 * coef[3] * 
                                           min(tm))
    }
    p0 <- exp(lp0)/(1 + exp(lp0))
    names(p0) <- NULL
    list(p0 = p0, extrap = extrap)
  }

# computes the functional anderson darling test
# f and g are two sets of functions, represented as matricies where each column is a separate function 
fadtest = function(f, g) {
  ### FAD method
  h = cbind(f, g)
  fpc = fpca.face(t(h))
  
  fscores = fpc$scores[1:ncol(f),]
  gscores = fpc$scores[(ncol(f) + 1):ncol(h),]
  
  fad = sapply(1:ncol(fscores), function(x) adk.test(fscores[,x], gscores[,x])$adk[2, 1:2])
  tvals = abs(fad[1,])
  pvals = fad[2,]
  
  c(max(tvals), min(1, min(pvals)*ncol(fscores)))
}

#### BAND test
# computes the BAND test 
# f and g are two sets of functions, represented as matricies where each column is a separate function 
bandtest = function(f, g) {
  ### Lopez and Romo
  fxd = MBD_relative(t(f), t(f))
  gxd = MBD_relative(t(g), t(f))
  
  q = mean(sapply(gxd, function(y) mean(fxd <= y)))
  pval = 1 - pnorm((0.5 - q) / sqrt(((1/ncol(f)) + (1/ncol(g)))/12))
  
  c(q, pval)
}
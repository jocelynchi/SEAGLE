#' Generate synthetic data according to a fixed effects model
#'
#' This function generates synthetic from the fixed effects model described in the
#' experimental studies portion of the paper.
#'
#' @param H Matrix of haplotype data (e.g. cosihap)
#' @param n Number of individuals
#' @param L Number of SNPs in the G matrix (Default is 100), should be a value between 1 and 604
#' @param maf Minor allele frequency (Default is 0.01)
#' @param gamma0 gamma0 Fixed effect coefficient for intercept (Default is 1)
#' @param gammaX gammaX Fixed effect coefficient for confounding covariates (Default is 1)
#' @param gammaE gammaE Fixed effect coefficient for E effect (Default is 1)
#' @param gammaG gammaG Fixed effect coefficient for G main effect
#' @param gammaGE gammaGE Fixed effect coefficient for GxE interaction effect
#' @param causal Number of causal SNPs (default is 40)
#' @param seed Seed (Default is 12345)
#'
#' @return Synthetic dataset containing y, X, E, G, epsilon, and number of causal SNPs
#'
#' @importFrom stats rnorm
#' @import Matrix
#'
#' @examples
#' dat <- makeSimData(H=cosihap, n=500, L=10, gammaG=1, gammaGE=0, causal=4, seed=1)
#'
#' @export
#'
makeSimData <- function(H, n, L=100, maf=0.01, gamma0=1, gammaX=1, gammaE=1, gammaG, gammaGE, causal=40, seed=12345) {

  N <- nrow(H)

  if (L <= 1 | L >= 604) {
    stop("L must be between 1 and 604.")
  }

  # Sample L of the MAF SNPs
  AF <- apply(H, 2, sum)/N
  MAF <- which(AF < maf)
  set.seed(seed)
  selects <- sample(MAF, L)
  Hnew <- H[, selects]

  # Make gammaGvec and gammaGEvec based on causal SNPs
  gammaGvec <- numeric(length=L)
  gammaGEvec <- numeric(length=L)
  gammaGvec[1:causal] <- gammaG
  gammaGEvec[1:causal] <- gammaGE

  # Make G
  G <- Matrix(0, nrow=n, ncol=L, sparse=TRUE)
  for (a in 1:n) {
    set.seed(seed+a)
    Hsample <- Hnew[sample(1:N, 2, replace=FALSE),]
    G[a,] <- apply(Hsample, 2, sum)
  }

  # Make y, X, E
  set.seed(seed*2)
  E <- as.numeric(scale(rnorm(n, 0, 1)))
  set.seed(seed*3)
  epsilon <- rnorm(n,0,1)
  set.seed(seed*4)
  X <- cbind(rep(1,n), scale(rnorm(n,0,1)))
  y <- X %*% c(gamma0, gammaX) + gammaE*E + G%*%gammaGvec + E*(G%*%gammaGEvec) + epsilon

  return(list(y=y, X=X, G=G, E=E, epsilon=epsilon, causal=causal))
}


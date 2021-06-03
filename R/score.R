#' Compute score-like test statistic and p-value for GxE test with SEAGLE algorithm
#'
#' This function computes the score test statistic and corresponding
#' p-value for the GxE test with the SEAGLE algorithm with input data
#' that have been prepared with the prep.SEAGLE function
#'
#' @param obj.SEAGLE Input data prepared with prep.SEAGLE function
#' @param init.tau Initial estimate for tau (Default is 0.5)
#' @param init.sigma Initial estimate for sigma (Default is 0.5)
#' @param pv Method of obtaining p-value (Either "liu" or "davies", Default is liu)
#'
#' @return Score-like test statistic T for the GxE effect and corresponding p-value
#'
#' @importFrom CompQuadForm davies liu
#'
#' @examples
#' dat <- makeSimData(H=cosihap, n=500, L=10, gammaG=1, gammaGE=0, causal=4, seed=1)
#' objSEAGLE <- prep.SEAGLE(y=dat$y, X=dat$X, intercept=1, E=dat$E, G=dat$G)
#' res <- SEAGLE(objSEAGLE, init.tau=0.5, init.sigma=0.5)
#'
#' @export
#'
SEAGLE <- function(obj.SEAGLE, init.tau=0.5, init.sigma=0.5, pv="liu") {

  y <- obj.SEAGLE$y
  Xtilde <- obj.SEAGLE$Xtilde
  qrXtilde <- obj.SEAGLE$qrXtilde
  E <- obj.SEAGLE$E
  G <- obj.SEAGLE$G

  vc <- estimate.vc(y, Xtilde, qrXtilde, beta, G, init.sigma, init.tau, tol=1e-3, maxiters=1000)
  tau <- vc$tau
  sigma <- vc$sigma

  n <- length(y)
  L <- dim(G)[2]

  # Pre-compute Vinv and Vhalf items
  tau_over_sigma <- tau/sigma
  GtG <- as.matrix(crossprod(G, G))   # LxL (small)
  M1 <- tau_over_sigma * GtG
  qrtausigGtG <- qr(M1)
  M <- diag(1, L) + M1
  qrM <- qr(M)

  # Compute test statistic
  XtVinvX <- t(Xtilde) %*% Vinv(G, qrM, tau_over_sigma, sigma, Xtilde)
  VinvY <- Vinv(G, qrM, tau_over_sigma, sigma, y)
  XtVinvY <- crossprod(Xtilde, VinvY)
  qrXtVinvX <- qr(XtVinvX)
  XtVinvXInvXtVinvY <- solve.qr(qrXtVinvX, XtVinvY)

  XT1 <- Xtilde %*% XtVinvXInvXtVinvY
  Py <- VinvY - Vinv(G, qrM, tau_over_sigma, sigma, XT1)

  Gtilde <- E * G
  TThalf <- crossprod(Gtilde, Py)
  teststat <- as.numeric((0.5) * crossprod(TThalf, TThalf))

  # Compute p-value
  VinvGtilde <- Vinv(G, qrM, tau_over_sigma, sigma, Gtilde)
  XtVinvGtilde <- crossprod(Xtilde, VinvGtilde)
  XT3 <- Xtilde %*% solve.qr(qrXtVinvX, XtVinvGtilde)
  T4 <- Vinv(G, qrM, tau_over_sigma, sigma, XT3)
  PGtilde <- VinvGtilde - T4
  C1tC1 <- t(Gtilde) %*% PGtilde / 2
  evC1tC1 <- eigen(C1tC1)$values
  evC <- evC1tC1[which(evC1tC1 > 1e-7)]

  if (pv == "davies") {
    pv.davies <- davies(teststat, evC[which(evC >= 1e-7)], sigma=sigma)$Qq
    return(list(T=teststat, pv=pv.davies))
  }
  else {
    pv.liu <- liu(teststat, evC[which(evC >= 1e-7)])
    return(list(T=teststat, pv=pv.liu))
  }
}

#' Function for applying V inverse in Algorithm 1
#'
#' This function applies V inverse via the Woodbury matrix identity
#'
#' @param G Matrix of genotype markers (size n x L)
#' @param qrM Pre-computation for LxL linear system solve
#' @param tau_over_sigma Tau over sigma from precomputation
#' @param sigma Variance component from model noise epsilon
#' @param RHS Matrix or vector on right-hand side of V inverse
#'
#' @return Matrix or vector resulting from left multiplication of Vinv with input RHS
#'
Vinv <- function(G, qrM, tau_over_sigma, sigma, RHS) {

  rhs = solve(qrM, as.matrix(t(G) %*% RHS))
  Out <- RHS - tau_over_sigma * (G %*% rhs)

  return((1/sigma)*Out)
}


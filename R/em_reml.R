#' REML EM Algorithm
#'
#' REML EM algorithm for estimating variance components
#'
#' @param y Vector of observed phenotypes
#' @param Xtilde Matrix of covariates (first column contains the intercept, last column contains the E factor for studying the GxE effect)
#' @param qrXtilde Object containing QR decomposition of Xtilde
#' @param beta Coefficient vector for covariate matrix Xtilde
#' @param G Matrix of genotype markers
#' @param init.sigma Initial sigma input (Default is 0.5)
#' @param init.tau Initial tau input (Default is 0.5)
#' @param tol Tolerance for convergence (Default is 1e-3)
#' @param maxiters Maximum number of iterations (Default is 1000)
#'
#' @return Estimates for tau and sigma
#'
estimate.vc <- function(y, Xtilde, qrXtilde, beta, G, init.sigma=0.5, init.tau=0.5, tol=1e-3, maxiters=1000) {

  # Prep data
  n <- length(y)
  L <- ncol(G)
  p <- qrXtilde$rank
  nminusp <- n - p
  tau_current <- init.tau
  sigma_current <- init.sigma

  # Prep fixed components
  u <- applyAt(qrXtilde, y)
  AtG <- applyAt(qrXtilde, G)
  GtA <- t(AtG)
  GtAu <- GtA %*% u
  GtAAtG <- GtA %*% AtG   # LxL (small)

  # Make containers for output
  sigmaDiff <- numeric(length=maxiters)
  tauDiff <- numeric(length=maxiters)

  for (a in 1:maxiters){

    # Compute sigma update
    temp1 <- sum((sigma_current * Rinv.u(G, AtG, GtAAtG, GtAu, u, tau_current, sigma_current))^2)
    temp2 <- sigma_current * tau_current * sum( diag(GtA %*% Rinv.AtG(G, AtG, GtAAtG, tau_current, sigma_current)))
    sigma_new <- (temp1 + temp2)/nminusp

    # Compute tau update
    temp1a <- GtA %*% Rinv.u(G, AtG, GtAAtG, GtAu, u, tau_current, sigma_current)
    temp1 <- tau_current^2 * sum(temp1a^2)
    temp2 <- tau_current^2 * GtA %*% Rinv.AtG(G, AtG, GtAAtG, tau_current, sigma_current)
    temp3 <- sum(tau_current - diag(temp2))
    tau_new <- (temp1 + temp3)/L

    # Check for convergence
    sigmaDiff[a] <- abs(sigma_current - sigma_new)/sigma_current
    tauDiff[a] <- abs(tau_current - tau_new)/tau_current
    if (sigmaDiff[a] < tol & tauDiff[a] < tol) { break }

    # Update sigma and tau
    sigma_current <- sigma_new
    tau_current <- tau_new
  }

  return(list(tau=tau_new, sigma=sigma_new))
}

#' Function for applying R inverse to u in REML EM algorithm
#'
#' Function for applying R inverse to u in REML EM algorithm
#'
#' @param G Matrix of genotype markers (size n x L)
#' @param AtG AtG from precomputation
#' @param GtAAtG GtAAtG from precomputation
#' @param GtAu GtAu from precomputation
#' @param u u=Aty from REML EM
#' @param tau Variance component from G main effect
#' @param sigma Variance component from model noise epsilon
#'
#' @return Vector resulting from left multiplication of Rinv with input vector u
#'
Rinv.u <- function(G, AtG, GtAAtG, GtAu, u, tau, sigma) {

  L <- dim(G)[2]
  tau_over_sigma <- tau / sigma       # needs to be updated in each EM iteration
  rhs = solve(diag(1, L) + tau_over_sigma * GtAAtG, GtAu)
  M <- u - tau_over_sigma * (AtG %*% rhs)

  return((1/sigma)*M)
}

#' Function for applying R inverse to AtG in REML EM algorithm
#'
#' Function for applying R inverse to AtG in REML EM algorithm
#'
#' @param G Matrix of genotype markers (size n x L)
#' @param AtG AtG from pre-computation
#' @param GtAAtG GtAAtG from pre-computation
#' @param tau Variance component from G main effect
#' @param sigma Variance component from model noise epsilon
#'
#' @return Matrix resulting from left multiplication of Rinv with input matrix AtG
#'
Rinv.AtG <- function(G, AtG, GtAAtG, tau, sigma) {

  L <- dim(G)[2]
  tau_over_sigma <- tau / sigma   # needs to be updated in each EM iteration
  rhs = solve(diag(1, L) + tau_over_sigma * GtAAtG, GtAAtG)
  M <- AtG - tau_over_sigma * (AtG %*% rhs)

  return((1/sigma)*M)
}


#' Function for applying t(A) on the left for REML EM
#'
#' @param qrXtilde Object from QR decomposition of Xtilde
#' @param RHS Object on right hand side of null of Xtilde^T
#'
#' @return Matrix or vector resulting from left multiplication of At with matrix or vector input RHS
#'
applyAt <- function(qrXtilde, RHS) {

  n <- dim(qrXtilde$qr)[1]
  p <- dim(qrXtilde$qr)[2]
  return(qr.qty(qrXtilde, as.matrix(RHS))[(p+1):n,])
}

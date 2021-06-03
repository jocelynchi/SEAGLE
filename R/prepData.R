#' Prepare data for input into SEAGLE function
#'
#' This function checks and formats data for input into SEAGLE function
#'
#' @param y Vector of observed phenotypes
#' @param X Matrix of covariates without genetic marker interactions
#' @param intercept 1 if the first column of X is the all ones vector, 0 otherwise
#' @param E E Vector of environmental covariates
#' @param G G Matrix of genotype data
#'
#' @return List object containing prepared data for input into SEAGLE function
#'
#' @examples
#' dat <- makeSimData(H=cosihap, n=500, L=10, gammaG=1, gammaGE=0, causal=4, seed=1)
#' objSEAGLE <- prep.SEAGLE(y=dat$y, X=dat$X, intercept=1, E=dat$E, G=dat$G)
#'
#' @importFrom stats lm
#' @export
#'
prep.SEAGLE <- function(y, X, intercept, E, G) {

  # Prep y
  y <- as.numeric(y)
  n <- length(y)

  # Prep X
  X <- as.matrix(X)
  if (dim(X)[1] != n) {
    stop("Please check your input data.  The phenotype vector and covariate matrix must have the same number of observations.")
  }

  # Prep E
  E <- as.numeric(E)
  if (length(E) != n) {
    stop("Please check your input data.  The environmental covariate must be a vector with the same number of rows as the phenotype vector.")
  }

  # Prep G
  G <- Matrix(G)
  if (dim(G)[1] != n) {
    stop("Please check your input data.  The phenotype vector and genetic marker matrix must have the same number of observations.")
  }

  # Prep Xtilde and qrXtilde
  if (intercept == 0) {
    Xtilde <- cbind(rep(1, n), X, scale(E))
  } else {
    Xtilde <- cbind(X, scale(E))
  }
  qrXtilde <- qr(Xtilde)

  # Get beta for REML EM
  beta <- lm(y ~ Xtilde - 1)

  # Return SEAGLE input object
  obj.SEAGLE <- list(y = y,
                     Xtilde = Xtilde,
                     E = E,
                     G = G,
                     beta = beta,
                     qrXtilde = qrXtilde
  )

  return(obj.SEAGLE)
}

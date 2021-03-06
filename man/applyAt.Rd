% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/em_reml.R
\name{applyAt}
\alias{applyAt}
\title{Function for applying t(A) on the left for REML EM}
\usage{
applyAt(qrXtilde, RHS)
}
\arguments{
\item{qrXtilde}{Object from QR decomposition of Xtilde}

\item{RHS}{Object on right hand side of null of Xtilde^T}
}
\value{
Matrix or vector resulting from left multiplication of At with matrix or vector input RHS
}
\description{
Function for applying t(A) on the left for REML EM
}

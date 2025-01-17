% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/est_MAP.R
\name{est.MAP}
\alias{est.MAP}
\title{Estimate traits based on genuine likelihood}
\usage{
est.MAP(
  FUN,
  responses,
  int,
  loads,
  uni,
  perms,
  SE = TRUE,
  lh.fun = lh,
  starts = NULL,
  box = 3,
  ...
)
}
\arguments{
\item{FUN}{function to compute response probability}

\item{responses}{matrix of block responses, rows = persons, columns = blocks. Responses should be given as indices for rank orders, corresponding to the columns in perms.}

\item{int}{vector of pair intercepts (i.e., intercepts for binary outcomes of pairwise comparisons)}

\item{loads}{matrix of item loadings, rows = items, columns = traits}

\item{uni}{matrix of item uniquenesses (diagonal)}

\item{perms}{matrix of permutations (i.e., rank orders). Can be obtained from calling permute()}

\item{SE}{logical. Obtain standard errors from generalized inverse of the negative hessian at the log-likelihood? defaults to TRUE.}

\item{lh.fun}{function to calculate likelihood across blocks. Defaults to lh.}

\item{starts}{matrix of starting values for the latent traits, rows = persons, columns = traits. If NULL, all starting values are zero. Defaults to NULL.}

\item{box}{numeric vector of length 1. Box constraints for the latent traits are set as $\pm$ box for all traits. Defaults to 3.}

\item{...}{additional arguments passed to FUN.}
}
\value{
list with 5 entries: traits = matrix of point estimates for the latent traits, row = persons, columns = traits. 
ses = matrix of standard errors for the trait estimates, if SE = FALSE, all entries are NA. 
errors, warns, messages = vectors of any errors, warnings and messages that occured during estimation, in the order of their occurence,
}
\description{
Estimate traits based on genuine likelihood
}

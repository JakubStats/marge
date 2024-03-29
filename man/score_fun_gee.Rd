% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MARGE_package.R
\name{score_fun_gee}
\alias{score_fun_gee}
\title{score_fun_gee}
\usage{
score_fun_gee(
  Y,
  N,
  n_vec,
  VS.est_list,
  AWA.est_list,
  J2_list,
  Sigma2_list,
  J11.inv,
  JSigma11,
  mu.est,
  V.est,
  B1,
  XA,
  nb = FALSE,
  ...
)
}
\arguments{
\item{Y}{: the response variable.}

\item{N}{: the number of clusters.}

\item{n_vec}{: a vector consisting of the cluster sizes for each cluster.}

\item{VS.est_list}{: a product of matrices.}

\item{AWA.est_list}{: a product of matrices.}

\item{J2_list}{: a product of matrices.}

\item{Sigma2_list}{: a product of matrices.}

\item{J11.inv}{: a product of matrices.}

\item{JSigma11}{: a product of matrices.}

\item{mu.est}{: estimates of the fitted mean under the null model.}

\item{V.est}{: estimates of the fitted variance under the null model.}

\item{B1}{: model matrix under the null model.}

\item{XA}{: model matrix under the alternative model.}

\item{nb}{: a logical argument, is the model a negative binomial model? The default is \code{FALSE}.}

\item{...}{: further arguments passed to or from other methods.}
}
\value{
\code{score_fun_gee} returns a calculated score statistic for the null and alternative model when fitting a GEE.
}
\description{
Given estimates from the null and the design matrix from alternative model, find the score statistic (this is used for GEEs only).
}
\references{
Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.

Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
}
\seealso{
\code{\link{score_fun_glm}}
}
\author{
Jakub Stoklosa and David I. Warton
}

#' stat_out_score_null
#'
#' A function that calculates parts for the score statistic for GEEs (it is used for the full path for forward selection).
#' @name stat_out_score_null
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param n : the maximum cluster size.
#' @param id : the ID for each individual in the cluster.
#' @param family : the specified "exponential" family for GLMs. The default is \code{family = "gaussian"}.
#' @param corstr : the specified "working correlation" structure. The default is \code{corstr = "independence"}.
#' @param B_null : model matrix under the null model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param is.gee : a logical argument, is this a GEE model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @details The null model used here is by definition the current parent model. We compare the alternative model (a new basis function combined with the parent) with the null model. In the code used, only the null model is fit (the dispersion parameter is also estimated for GEE) and then the score statistic is obtained.
#' @return \code{stat_out_score_null} returns a list of values (mainly products of matrices) that make up the final score statistic calculation (required for another function).
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @importFrom stats binomial poisson
#' @export
#' @seealso \code{\link{stat_out}} and \code{\link{stat_out_score_glm_null}}
stat_out_score_null <- function(Y, N, n, id, family = "gaussian", corstr = "independence", B_null, nb = FALSE, is.gee = FALSE, ...) {
  n_vec <- as.numeric(table(id))

  if (is.gee) {
    # if (family == "gaussian") {
    #   ests <- geepack::geeglm(Y ~ B_null - 1, id = id, corstr = corstr, ...)
    #   alpha.est <- ests$geese$alpha
    # }
    # if (family != "gaussian") {
    #   ests <- geepack::geeglm(Y ~ B_null - 1, id = id, family = family, corstr = corstr, ...)
    #   alpha.est <- ests$geese$alpha
    # }
    # ED: we can just pass the family argument and avoid the above switch
    ests <- geepack::geeglm(Y ~ B_null - 1, id = id, family = family, corstr = corstr, ...)
    alpha.est <- ests$geese$alpha

    mu.est <- as.matrix(stats::fitted.values(ests))
  }

  if (!is.gee) {
    if (nb == TRUE) {
      ests <- gamlss::gamlss(Y ~ B_null - 1, family = "NBI", trace = FALSE)
      mu.est <- as.matrix(stats::fitted.values(ests))
    }
    if (nb == FALSE) {
      # if (family == "gaussian") ests <- stats::glm.fit(B_null, Y, ...)
      # if (family == "binomial") ests <- stats::glm.fit(B_null, Y, family = binomial(link = "logit"), ...)
      # if (family == "poisson") ests <- stats::glm.fit(B_null, Y, family = poisson(link = "log"), ...)
      # ED: we can just pass the family argument (if we move to glm() and note that these use only the default links anyway) and avoid the above switch
      ests <-stats::glm(Y ~ B_null - 1, family = family)
      mu.est <- as.matrix(stats::fitted.values(ests))
    }
    alpha.est <- 1
  }

  # if (family == "gaussian") V.est <- rep(1, nrow(mu.est))
  # if (family == "binomial") V.est <- mu.est*(1 - mu.est)
  # if (family == "poisson") {
  #   if (nb == FALSE) V.est <- mu.est
  #   if (nb == TRUE) V.est <- mu.est*(1 + mu.est*(exp(ests$sigma.coef)))
  # }
  # ED: a switch function may be more efficient here
  V.est <- switch(family,
                  gaussian = rep(1, nrow(mu.est)),
                  binomial = mu.est*(1 - mu.est),
                  poisson = if (nb) {mu.est*(1 + mu.est*(exp(ests$sigma.coef)))} else {mu.est}
                  )

  p <- ncol(B_null)
  n_vec1 <- c(0, n_vec)

  VS.est_list <- list()
  AWA.est_list <- list()
  J2_list <- list()
  Sigma2_list <- list()

  J11 <- matrix(0, nrow = p, ncol = p)
  Sigma11 <- matrix(0, nrow = p, ncol = p)

  for (i in 1:N) {
    k <- sum(n_vec[1:i])
    # if (corstr == "independence") R_alpha <- diag(1, nrow = n_vec[i], ncol = n_vec[i])
    # if (corstr == "exchangeable") R_alpha <- matrix(c(rep(alpha.est, n_vec[i]*n_vec[i])), ncol = n_vec[i]) + diag(c(1 - alpha.est), ncol = n_vec[i], nrow = n_vec[i])
    # if (corstr == "ar1") R_alpha <- alpha.est^outer(1:n_vec[i], 1:n_vec[i], function(x, y) abs(x - y))
    # ED: again swapping this for a switch function
    R_alpha <- switch(corstr,
                      independence = diag(1, nrow = n_vec[i], ncol = n_vec[i]),
                      exchangeable = matrix(c(rep(alpha.est, n_vec[i]*n_vec[i])), ncol = n_vec[i]) + diag(c(1 - alpha.est), ncol = n_vec[i], nrow = n_vec[i]),
                      ar1 = alpha.est^outer(1:n_vec[i], 1:n_vec[i], function(x, y) abs(x - y))
                      )

    V.est_i <- diag(sqrt(V.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])%*%R_alpha%*%diag(sqrt(V.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])
    V.est_i_inv <- chol2inv(chol(V.est_i))
    # ED: these can be a bottleneck though speeding them up requires some extra thought/time (note to self: Simon Wood book for tips)

    # S.est_i <- c(t(Y))[(sum(n_vec1[1:i]) + 1):k] - mu.est[(sum(n_vec1[1:i]) + 1):k]
    # ED: can remove the transposing of Y
    S.est_i <- Y[(sum(n_vec1[1:i]) + 1):k] - mu.est[(sum(n_vec1[1:i]) + 1):k]

    # AWA.est_i <- V.est_i_inv%*%(S.est_i%*%t(S.est_i))%*%V.est_i_inv
    # ED: removing the transpose by going straight to the outer product
    AWA.est_i <- V.est_i_inv%*%(S.est_i%o%S.est_i)%*%V.est_i_inv

    # if (nb == FALSE) D.est_i <- diag((V.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])%*%B_null[(sum(n_vec1[1:i]) + 1):k, ]
    # if (nb == TRUE) D.est_i <- diag((mu.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])%*%B_null[(sum(n_vec1[1:i]) + 1):k, ]
    # ED: removing double if check
    if (nb) {
      D.est_i <- diag((mu.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])%*%B_null[(sum(n_vec1[1:i]) + 1):k, ]
    } else {
      D.est_i <- diag((V.est[(sum(n_vec1[1:i]) + 1):k]), nrow = n_vec[i], ncol = n_vec[i])%*%B_null[(sum(n_vec1[1:i]) + 1):k, ]
    }

    # ED: want to transpose D.est_i just once (as it is used a few times ahead)
    tD.est_i <- t(D.est_i)

    # J1_i <- t(D.est_i)%*%V.est_i_inv%*%D.est_i
    # J11 <- J11 + J1_i
    # J2_i <- t(D.est_i)%*%V.est_i_inv
    # ED: can just calc J2_i once and use it for J1_i (also can we remove the transpose by flipping the order?)
    J2_i <-  tD.est_i%*%V.est_i_inv
    J1_i <- J2_i%*%D.est_i
    J11 <- J11 + J1_i

    # Sigma1_i <- t(D.est_i)%*%AWA.est_i%*%(D.est_i)
    # Sigma11 <- Sigma11 + Sigma1_i
    # Sigma2_i <- t(D.est_i)%*%AWA.est_i
    # ED: again can just calc Sigma2_i once and use it for Sigma1_i
    Sigma2_i <- tD.est_i%*%AWA.est_i
    Sigma1_i <- Sigma2_i%*%D.est_i
    Sigma11 <- Sigma11 + Sigma1_i

    VS.est_list <- c(VS.est_list, list(V.est_i_inv%*%S.est_i))
    AWA.est_list <- c(AWA.est_list, list(AWA.est_i))

    J2_list <- c(J2_list, list(J2_i))
    Sigma2_list <- c(Sigma2_list, list(Sigma2_i))
  }

  J11.inv <- chol2inv(chol(J11))

  JSigma11 <- J11.inv%*%Sigma11%*%J11.inv

  list(VS.est_list = VS.est_list, AWA.est_list = AWA.est_list, J2_list = J2_list, J11.inv = J11.inv, Sigma2_list = Sigma2_list, JSigma11 = JSigma11, mu.est = mu.est, V.est = V.est)
}

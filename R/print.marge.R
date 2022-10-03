#' Print for objects of class 'marge'
#'
#' @param x a marge model object
#' @param pen a set penalty used for the GCV (which model information to return). The default is "2".
#' @param ...
#'
#' @return printed argument
#' @exportS3Method base::print marge
#'
#' @examples
#' #' data(leptrine)
#' dat1 <- leptrine[[1]]
#' mod <- marge(Y ~ RAIN_DRY_QTR + FC, data = dat1, family = "binomial", N = nrow(dat1), n = 1)
#' mod
print.marge <- function(x, ...) {

  mod <- x$final_mod
  # get the coefficients
  tmp.coefs <- data.frame(mod$coefficients)
  dimnames(tmp.coefs)[[2]] <- dimnames(attr(x$terms, "factors"))[[1]][1]
  # get the summaries
  # get the summaries
  if (x$is.gee) {
    tmp.stats <- with(mod, c(corstr, table(geese$clusz), round(geese$gamma, 3)))
    names(tmp.stats) <- c("Correlation Structure", "Clusters", "Estimated Scale Parameters:")
  } else {
    tmp.stats <- with(mod, c(round(null.deviance, 2), round(df.null, 0), round(deviance, 2), round(df.residual, 0), round(deviance/null.deviance, 2), round(aic, 2), round(iter, 0), round(converged, 0)))
    names(tmp.stats) <- c("nulldev", "df", "dev", "df", "devratio", "AIC", "iters", "converged")
  }

  # print the output
  cat("Call: ", deparse(x$call), "\n\nCoefficients:\n")
  print(mod$family)
  print(tmp.stats)
  cat("MARGE GCV ", x$GCV)
  # cat("Null dev: ", mod$deviance, " Residual deviance: " mod)

}

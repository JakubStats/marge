#' marge
#'
#' MARS fitting function for generalized linear models (GLM) and generalized estimating equations (GEE).
#' @name marge
#' @param formula : an object of class "formula" (or one that can be coerced to that class): a symbolic description of first order terms to be included in the model to be fitted. See GLM function for further formula details.
#' @param data : a data frame containing all the terms found in \code{formula}. Should have n by N rows.
#' @param N : the number of clusters.
#' @param n : the maximum cluster size.
#' @param id : a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
#' @param family : the specified family for the GLM/GEE. The default is \code{family = "gaussian"}. The current available families are: "gaussian", "binomial", "poisson" and the "negative binomial" (to get that use family = "poisson" and nb = TRUE).
#' @param corstr : the specified "working correlation" structure for the GEE. The default is \code{corstr = "independence"}.
#' @param pen : a character string describing the penalty used for the GCV (note: MARGE doesn't actually use this). The default is 2.
#' @param tols_score : the set tolerance for monitoring the convergence for the difference in score statistics between the parent and candidate model (this is the lack-of-fit criterion used for MARGE). The default is 0.00001
#' @param M : a set threshold for the number of basis functions to be used. The default is 21.
#' @param minspan : a set minimum span value. The default is \code{minspan = NULL}.
#' @param print.disp : a logical argument, do you want to print the output? The default is \code{FALSE}.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param is.gee : is this a GEE model? The default is \code{FALSE}.
#' @param plot.wic : a logical indicating whether the WIC should be plotted for the backward path. The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @details For further details please look at the \code{mars_ls} function - there are more details on the general MARS algorithm. MARGE will produce output for two penalties: 2 and log(N). A figure is automatically generated plotting WIC against the no. of parameters.
#' @return \code{marge} returns a list of calculated values consisting of:
#' @return \code{B_final} : the basis model matrix for the final model fit.
#' @return \code{wic_mat} : a matrix of WIC values (with both penalties) for MARGE models given by the forward pass.
#' @return \code{min_wic} : the WIC (with both penalties) for the final MARGE model.
#' @return \code{GCV} : the GCV for the final selected model.
#' @return \code{y_pred} : the fitted values from the final selected model (with both penalties).
#' @return \code{final_mod} : the final selected (with both penalties) model matrix.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, \strong{19}, 1--67.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, \strong{70}, 110--120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @noRd
#' @seealso \code{\link{mars_ls}} and \code{\link{backward_sel_WIC}}
#' @importFrom gamlss gamlss
#' @importFrom geepack geeglm
#' @importFrom mvabund manyglm
#' @importFrom MASS glm.nb
#' @importFrom stats binomial poisson
#' @examples
#' # Example 1: Tree growth data from Wood (2006) monograph
#'
#' library(mgcv)
#'
#' N <- 14                            # Sample size (number of clusters).
#' n <- 6                             # Cluster size.
#' id <- factor(rep(1:N, each = n))   # The ID of each cluster.
#'
#' Y0 <- matrix(log(Loblolly$height), ncol = n, byrow = TRUE) # Response variable.
#' Y <- c(t(Y0))
#' X_pred0 <- matrix(log(Loblolly$age), ncol = n, byrow = TRUE) # Predictor variable.
#' X_pred <- scale(c(t(X_pred0)))  # Standarize the predictor.
#' colnames(X_pred) <- c("age")
#'
#' dat1 <- as.data.frame(cbind(X_pred, c(t(Y)), id))
#' colnames(dat1) <- c("age", "Y", "id")
#' dat1$age1 <- c(t(X_pred0))
#' dat1$id <- factor(dat1$id)
#'
#' family <- "gaussian"   # The selected "exponential" family for the GLM/GEE.
#'
#' is.gee <- TRUE         # Is the model a GEE?
#' nb <- FALSE            # Is this a negative binomial model?
#'
#' tols_score <- 0.00001  # Tolerance for stopping condition in forward pass for GEE.
#' M <- 21                # Max. no. of terms.
#' pen <- 2               # Penalty to be used in GCV.
#' minspan <- NULL        # A set minimum span value.
#' print.disp <- FALSE    # Print ALL the output?
#'
#' # Start model fitting functions here.
#'
#' corstr <- "independence"  # Independent working correlation structure.
#' model_marge_ind <- marge(X_pred, Y, N, n, id, family, corstr, pen, tols_score,
#' M, minspan, print.disp, nb, is.gee)
#'
#' corstr <- "ar1"           # AR1 working correlation structure.
#' model_marge_ar1 <- marge(X_pred, Y, N, n, id, family, corstr, pen, tols_score,
#' M, minspan, print.disp, nb, is.gee)
#'
#' corstr <- "exchangeable"  # Exchangeable working correlation structure.
#' model_marge_exch <- marge(X_pred, Y, N, n, id, family, corstr, pen, tols_score,
#' M, minspan, print.disp, nb, is.gee)
#'
#' # Example 2: Presence-absence data
#'
#' # Load the "leptrine" presence-absence data.
#'
#' data(leptrine)
#'
#' dat1 <- leptrine[[1]]  # Training data.
#' Y <- dat1$Y            # Response variable.
#' N <- length(Y)         # Sample size (number of clusters).
#' n <- 1                 # Cluster size.
#' id <- rep(1:N, each = n)  # The ID of each cluster.
#'
#' X_pred <- dat1[, -c(3:10)] # Design matrix using only two (of nine) predictors.
#'
#' # Set MARGE tuning parameters.
#'
#' family <- "binomial"   # The selected "exponential" family for the GLM/GEE.
#' is.gee <- FALSE        # Is the model a GEE?
#' nb <- FALSE            # Is this a negative binomial model?
#'
#' tols_score <- 0.0001   # A set tolerance (stopping condition) in forward pass for MARGE.
#' M <- 21                # A set threshold for the maximum number of basis functions to be used.
#' print.disp <- FALSE    # Print ALL the output?
#' pen <- 2               # Penalty to be used in GCV.
#' minspan <- NULL        # A set minimum span value.
#'
#' # Fit the MARGE models (about ~ 30 secs.)
#'
#' mod <- marge(Y ~ RAIN_DRY_QTR + FC, data = dat1, N, n, id, family, corstr, pen, tols_score,
#'              M, minspan, print.disp, nb, is.gee)
marge2 <- function (formula = formula(data), sformula = ~1, id, waves = NULL,
                   data = parent.frame(), subset = NULL, na.action = na.omit,
                   contrasts = NULL, weights = NULL, zcor = NULL, corp = NULL,
                   control = geese.control(...), b = NULL, alpha = NULL, gm = NULL,
                   family = gaussian(), mean.link = NULL, variance = NULL, cor.link = "identity",
                   sca.link = "identity", link.same = TRUE, scale.fix = FALSE,
                   scale.value = 1, corstr = "independence", ...)
{
  scall <- match.call()
  mnames <- c("", "formula", "data", "offset", "weights", "subset",
              "na.action", "id", "waves", "corp")
  cnames <- names(scall)
  cnames <- cnames[match(mnames, cnames, 0)]
  mcall <- scall[cnames]
  if (is.null(mcall$id))
    mcall$id <- as.name("id")
  mcall[[1]] <- as.name("model.frame")
  m <- eval(mcall, parent.frame())
  y <- model.extract(m, "response")
  if (is.null(dim(y)))
    N <- length(y)
  else N <- dim(y)[1]
  mterms <- attr(m, "terms")
  x <- model.matrix(mterms, m, contrasts)
  offset <- model.extract(m, "offset")
  if (is.null(offset))
    offset <- rep(0, N)
  w <- model.extract(m, "weights")
  if (is.null(w))
    w <- rep(1, N)
  id <- model.extract(m, id)
  waves <- model.extract(m, "waves")
  corp <- model.extract(m, "corp")
  if (is.null(id))
    stop("id variable not found.")
  mcall$formula <- formula
  mcall$formula[3] <- switch(match(length(sformula), c(0, 2, 3)), 1, sformula[2], sformula[3])
  m <- eval(mcall, parent.frame())
  terms <- attr(m, "terms")
  zsca <- model.matrix(terms, m, contrasts)
  soffset <- model.extract(m, "offset")
  if (is.null(soffset))
    soffset <- rep(0, N)
  if (is.character(family))
    family <- get(family)
  if (is.function(family))
    family <- family()
  ans <- list(x, y, id, offset, soffset, w, waves, zsca,
                   zcor, corp, control, b, alpha, gm, family, mean.link,
                   variance, cor.link, sca.link, link.same, scale.fix, scale.value,
                   corstr, ...)
  ans <- c(ans, list(call = scall, formula = formula))
  # class(ans) <- "geese"
  m
}

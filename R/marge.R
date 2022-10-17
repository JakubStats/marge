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
#' @export
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
marge <- function(formula, data, N, n = 1, id = c(1:nrow(data)), family = "gaussian", corstr = "independence", pen = c("two", "log"), tols_score = 0.00001, M = 21, minspan = NULL, print.disp = FALSE, nb = FALSE, is.gee = FALSE, plot.wic = FALSE, ...) {

  # changing 'pen' argument
  pen <- match.arg(pen)
  preds <- all.vars(formula[[3]])
  resp <- all.vars(formula[[2]])
  if (!all(preds %in% colnames(data))) {
    stop("Not all fixed effect terms found in the data provided")
  }
  if (!resp %in% colnames(data)) {
    stop("Response term not found in the data provided")
  }
  if (!is.logical(is.gee)) {
    stop("'is.gee' must be either TRUE or FALSE")
  }
  X_pred <- data.frame(data[ , preds])
  Y <- as.vector(data[ , resp])

  NN <- length(Y)    # Total sample size = N*n.

  q <- ncol(X_pred)  # No. of predictor variables.

  n_vec <- as.numeric(table(id))

  # Algorithm 2 (forward pass) as in Friedman (1991). Uses score statistics instead of RSS, etc.

  B <- as.matrix(rep(1, NN))       # Start with the intercept model.

  colnames(B) <- c("Intercept")
  var_name_vec <- c("Intercept")
  var_name_list <- list("Intercept")
  B_names_vec <- c("Intercept")
  min_knot_vec <- c("Intercept")
  pred.name_vec <- c("Intercept")
  cut_vec <- c("Intercept")
  trunc.type_vec <- c(1)
  is.int_vec <- c("FALSE")
  mod_struct <- c(1)          # Univariate (1) or interaction (2).

  score_term <- c(0)

  TSS <- sum((Y - mean(Y))^2)
  GCV.null <- TSS/(NN*(1 - 1/NN)^2)

  # Null model setup.

  m <- 1
  k <- 1
  breakFlag <- FALSE
  ok <- TRUE
  int.count <- 0

  while(ok) {
    if (breakFlag) break

    var.mod_temp <- c()
    score_term_temp <- c()
    min_knot_vec_temp <- c()
    int.count1_temp <- c()
    is.int_temp <- c()
    trunc.type_temp <- c()

    B_new_list_temp <- list()
    var_name_list1_temp <- list()
    B_names_temp <- list()
    X_red_temp <- list()
    B_temp_list <- list()

    # Obtain/calculate the null stats here (speeds things up).

    if (is.gee | n > 1) {
      B_null_stats <- stat_out_score_null(Y, N, n, id, family, corstr, B, nb, is.gee, ...) # ED: I have done some optimisation of this

      VS.est_list <- B_null_stats$VS.est_list
      AWA.est_list <- B_null_stats$AWA.est_list
      J2_list <- B_null_stats$J2_list
      J11.inv <- B_null_stats$J11.inv
      Sigma2_list <- B_null_stats$Sigma2_list
      JSigma11 <- B_null_stats$JSigma11
      mu.est <- B_null_stats$mu.est
      V.est <- B_null_stats$V.est
    }
    if (!is.gee & n == 1) {
      B_null_stats <- stat_out_score_glm_null(Y, family, B, nb, ...) # ED: seems quicker in this setting though

      VS.est_list <- B_null_stats$VS.est_list
      A_list <- B_null_stats$A_list
      B1_list <- B_null_stats$B1_list
      mu.est <- B_null_stats$mu.est
      V.est <- B_null_stats$V.est
    }

    for (v in 1:q) {
      var_name <- colnames(X_pred)[v]
      X <- round(X_pred[, v], 4)

      X_red1 <- min_span(X, q, minspan)  # Reduce the space between knots.
      X_red2 <- max_span(X, q)  # Truncate the ends of data to avoid extreme values.

      X_red <- intersect(X_red1, X_red2)

      score_knot_both_int_mat <- c()
      score_knot_both_add_mat <- c()
      score_knot_one_int_mat <- c()
      score_knot_one_add_mat <- c()

      int.count1 <- 0

      if (ncol(B) > 1) {
        in.set <- sum(!var_name_vec%in%var_name)
      } else {
        in.set <- 0
      }

      for (t in 1:length(X_red)) {
        b1_new <- matrix(tp1(X, X_red[t]), ncol = 1)  # Pairs of truncated functions.
        b2_new <- matrix(tp2(X, X_red[t]), ncol = 1)

        score_knot_both_int <- c()
        score_knot_both_add <- c()
        score_knot_one_int <- c()
        score_knot_one_add <- c()

        # ED: work out the "in.set > 0" scenario first as aspects of "in.set == 0" are shared by both
        if (in.set > 0) {
          var_name_struct <- which(((var_name != var_name_vec)*mod_struct) == 1)
          colnames(B)[1] <- c("")
          B2 <- as.matrix(B[, var_name_struct])
          if (k != 1 & (sum(!var_name_vec[-1]%in%var_name) > 0)) B2 <- as.matrix(B2[, -1])
          if (ncol(B2) == 0) B2 <- as.matrix(B[, 1])

          for (nn in 1:ncol(B2)) {
            B2a <- matrix(rep(B2[, nn], 2), ncol = 2)
            B2b <- matrix(B2[, nn], ncol = 1)
            B_new_both_int <- cbind(B, B2a*cbind(b1_new, b2_new))
            B_new_one_int <- cbind(B, B2b*b1_new)     # Interaction model with one truncated function (i.e., the positive part).

            if (is.gee == TRUE | n > 1) meas_model_both_int <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_both_int, B2a*cbind(b1_new, b2_new), nb, ...)
            if (is.gee == FALSE & n == 1) meas_model_both_int <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_both_int, B2a*cbind(b1_new, b2_new), nb, ...)
            if (is.gee == TRUE | n > 1) meas_model_one_int <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_one_int, B2b*b1_new, nb, ...)
            if (is.gee == FALSE & n == 1) meas_model_one_int <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_one_int, B2b*b1_new, nb, ...)

            score_knot_both_int <- c(score_knot_both_int, meas_model_both_int$score)
            score_knot_one_int <- c(score_knot_one_int, meas_model_one_int$score)
          }
        } else  if (in.set == 0) {
          score_knot_both_int = score_knot_one_int = -100000   # Interaction set is impossible since there is nothing to interact with, so let the LOF measure be a huge negative number.
        }

        B_new_both_add <- cbind(B, b1_new, b2_new)   # Additive model with both truncated functions.
        B_new_one_add <- cbind(B, b1_new)       # Additive model with one truncated function (positive part).

        # ED: reducing # if statements
        if (is.gee == TRUE | n > 1) {
          meas_model_both_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
          meas_model_one_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        } else if (is.gee == FALSE & n == 1) {
          meas_model_both_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
          meas_model_one_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        }

        score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
        score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)
        # ED: this last part of the code is the same for either scenario - compare to commented out code below

        # if (in.set == 0) {
        #   B_new_both_add <- cbind(B, b1_new, b2_new)   # Additive model with both truncated functions.
        #   B_new_one_add <- cbind(B, b1_new)       # Additive model with one truncated function (positive part).
        #
        #   if (is.gee == TRUE | n > 1) meas_model_both_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
        #   if (is.gee == FALSE & n == 1) meas_model_both_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
        #   if (is.gee == TRUE | n > 1) meas_model_one_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        #   if (is.gee == FALSE & n == 1) meas_model_one_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        #
        #   score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
        #   score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)
        #
        #   score_knot_both_int = score_knot_one_int = -100000   # Interaction set is impossible since there is nothing to interact with, so let the LOF measure be a huge negative number.
        # }
        #
        # if (in.set > 0) {
        #   var_name_struct <- which(((var_name != var_name_vec)*mod_struct) == 1)
        #   colnames(B)[1] <- c("")
        #   B2 <- as.matrix(B[, var_name_struct])
        #   if (k != 1 & (sum(!var_name_vec[-1]%in%var_name) > 0)) B2 <- as.matrix(B2[, -1])
        #   if (ncol(B2) == 0) B2 <- as.matrix(B[, 1])
        #
        #   for (nn in 1:ncol(B2)) {
        #     B2a <- matrix(rep(B2[, nn], 2), ncol = 2)
        #     B2b <- matrix(B2[, nn], ncol = 1)
        #     B_new_both_int <- cbind(B, B2a*cbind(b1_new, b2_new))
        #     B_new_one_int <- cbind(B, B2b*b1_new)     # Interaction model with one truncated function (i.e., the positive part).
        #
        #     if (is.gee == TRUE | n > 1) meas_model_both_int <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_both_int, B2a*cbind(b1_new, b2_new), nb, ...)
        #     if (is.gee == FALSE & n == 1) meas_model_both_int <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_both_int, B2a*cbind(b1_new, b2_new), nb, ...)
        #     if (is.gee == TRUE | n > 1) meas_model_one_int <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_one_int, B2b*b1_new, nb, ...)
        #     if (is.gee == FALSE & n == 1) meas_model_one_int <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_one_int, B2b*b1_new, nb, ...)
        #
        #     score_knot_both_int <- c(score_knot_both_int, meas_model_both_int$score)
        #     score_knot_one_int <- c(score_knot_one_int, meas_model_one_int$score)
        #   }
        #
        #   B_new_both_add <- cbind(B, b1_new, b2_new)
        #   B_new_one_add <- cbind(B, b1_new)
        #
        #   if (is.gee == TRUE | n > 1) meas_model_both_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
        #   if (is.gee == FALSE & n == 1) meas_model_both_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_both_add, cbind(b1_new, b2_new), nb, ...)
        #   if (is.gee == TRUE | n > 1) meas_model_one_add <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        #   if (is.gee == FALSE & n == 1) meas_model_one_add <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_new_one_add, b1_new, nb, ...)
        #
        #   score_knot_both_add <- c(score_knot_both_add, meas_model_both_add$score)
        #   score_knot_one_add <- c(score_knot_one_add, meas_model_one_add$score)
        # }

        score_knot_both_int_mat <- rbind(score_knot_both_int_mat, score_knot_both_int)
        score_knot_both_add_mat <- rbind(score_knot_both_add_mat, score_knot_both_add)
        score_knot_one_int_mat <- rbind(score_knot_one_int_mat, score_knot_one_int)
        score_knot_one_add_mat <- rbind(score_knot_one_add_mat, score_knot_one_add)
      }

      # See the LM code above in regards to what the conditions below do.

      if (sum(!(apply(score_knot_both_int_mat, 1, is.na))) == 0 & sum(!(apply(score_knot_one_int_mat, 1, is.na))) == 0) {
        int <- FALSE
        if (sum(!is.na(score_knot_both_add_mat)) > 0 & sum(!is.na(score_knot_one_add_mat)) > 0) {
          if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 2
            score_knot <- score_knot_both_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
          }
          if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            score_knot <- score_knot_one_add_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
          }
        }
        if (sum(!is.na(score_knot_both_add_mat)) == 0 & sum(!is.na(score_knot_one_add_mat)) > 0) {
          trunc.type <- 1
          score_knot <- score_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
        }
        if (sum(!is.na(score_knot_both_add_mat)) > 0 & sum(!is.na(score_knot_one_add_mat)) == 0) {
          trunc.type <- 2
          score_knot <- score_knot_one_add_mat
          min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
        }
        if (sum(!is.na(score_knot_both_add_mat)) == 0 & sum(!is.na(score_knot_one_add_mat)) == 0) {
          breakFlag <- TRUE
          break
        }
      }

      if (sum(!(apply(score_knot_both_int_mat, 1, is.na))) == 0 & sum(!(apply(score_knot_one_int_mat, 1, is.na))) > 0) {
        if (sum(!is.na(score_knot_both_add_mat)) == 0) {
          trunc.type <- 1
          if (sum(!is.na(score_knot_one_add_mat)) == 0) {
            int <- TRUE
            score_knot <- score_knot_one_int_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
            best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
          }
          if (sum(!is.na(score_knot_one_add_mat)) > 0) {
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        }
        if (sum(!is.na(score_knot_both_add_mat)) > 0) {
          if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 1
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
            if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
          }
          if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE)) <= utils::tail(max(score_knot_both_add_mat, na.rm = TRUE))) {
            int <- FALSE
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 2
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        }
      }

      if (sum(!(apply(score_knot_both_int_mat, 1, is.na))) > 0 & sum(!(apply(score_knot_one_int_mat, 1, is.na))) == 0) {
        if (sum(!is.na(score_knot_both_add_mat)) == 0) {
          if (sum(!is.na(score_knot_one_add_mat)) == 0) {
            int <- TRUE
            trunc.type <- 2
            score_knot <- score_knot_both_int_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
            best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
          }
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- TRUE
            trunc.type <- 2
            score_knot <- score_knot_both_int_mat
            min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
            best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
          }
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            trunc.type <- 1
            score_knot <- score_knot_one_add_mat
            min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
          }
        }
        if (sum(!is.na(score_knot_both_add_mat)) > 0) {
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            trunc.type <- 2
            if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              score_knot <- score_knot_both_int_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
            if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
            }
          }
          if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1)) {
            int <- FALSE
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 2
              score_knot <- score_knot_both_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
          }
        }
      }

      if (sum(!(apply(score_knot_both_int_mat, 1, is.na))) > 0 & sum(!(apply(score_knot_one_int_mat, 1, is.na))) > 0) {
        if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
          if (sum(!is.na(score_knot_both_add_mat)) > 0) {
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 2
                score_knot <- score_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)))
              }
            }
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                trunc.type <- 2
                score_knot <- score_knot_both_int_mat
                min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
                best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
              }
              if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            }
          }
          if (sum(!is.na(score_knot_both_add_mat)) == 0) {
            if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              trunc.type <- 1
              score_knot <- score_knot_one_add_mat
              min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
            }
            if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1)) {
              int <- TRUE
              trunc.type <- 2
              score_knot <- score_knot_both_int_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
          }
        }
        if (utils::tail(max(score_knot_both_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
          if (sum(!is.na(score_knot_both_add_mat)) > 0) {
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
              int <- FALSE
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 1
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
              if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                trunc.type <- 2
                score_knot <- score_knot_both_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            }
            if (utils::tail(max(score_knot_both_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
              trunc.type <- 1
              if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) > utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                score_knot <- score_knot_one_int_mat
                min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
                best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
              }
              if (utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1) <= utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)), n = 1)
              }
            }
          }
          if (sum(!is.na(score_knot_both_add_mat)) == 0) {
            trunc.type <- 1
            if (sum(!is.na(score_knot_one_add_mat)) == 0) {
              int <- TRUE
              score_knot <- score_knot_one_int_mat
              min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
              best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
            }
            if (sum(!is.na(score_knot_one_add_mat)) > 0) {
              if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) >= utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
                int <- FALSE
                score_knot <- score_knot_one_add_mat
                min_knot1 <- utils::tail(which.max(round(score_knot, 6)))
              }
              if (utils::tail(max(score_knot_one_add_mat, na.rm = TRUE), n = 1) < utils::tail(max(score_knot_one_int_mat, na.rm = TRUE), n = 1)) {
                int <- TRUE
                score_knot <- score_knot_one_int_mat
                min_knot1 <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[1]
                best.var <- utils::tail(which(utils::tail(max(round(score_knot, 6), na.rm = TRUE), n = 1) == round(score_knot, 6), arr.ind = TRUE), n = 1)[2]
              }
            }
          }
        }
      }

      b1_new <- matrix(tp1(X, X_red[min_knot1]), ncol = 1)
      b2_new <- matrix(tp2(X, X_red[min_knot1]), ncol = 1)
      colnames(b1_new) <- var_name
      colnames(b2_new) <- var_name

      B_name1 <- paste("(", var_name, "-", signif (X_red[min_knot1], 4), ")", sep = "")
      B_name2 <- paste("(", signif (X_red[min_knot1], 4), "-", var_name, ")", sep = "")

      if (int == TRUE) {
        mod_struct1 <- which(mod_struct == 1)
        colnames(B)[1] <- c("")
        var_name1 <- which(var_name_vec != var_name)
        if (int.count == 0) var_name2 <- var_name_vec[var_name1]
        if (int.count > 0) var_name2 <- var_name_vec[var_name_struct]
        var_name_struct <- mod_struct1[mod_struct1%in%var_name1]
        B2 <- as.matrix(B[, var_name_struct])
        B3_names <- B_names_vec[var_name_struct]
        B3_names <- B3_names[-1]
        B2 <- as.matrix(B2[, -1])
        var_name2 <- var_name2[-1]

        if (trunc.type == 2) {
          B2a <- matrix(rep(B2[, best.var], 2), ncol = 2)
          B_temp <- cbind(B, B2a*cbind(b1_new, b2_new))
          B_new <- B2a*cbind(b1_new, b2_new)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 2)
        }

        if (trunc.type == 1) {
          B2b <- matrix(B2[, best.var], ncol = 1)
          B_temp <- cbind(B, B2b*b1_new)     # Interaction model with one truncated basis function (i.e., the positive part).
          B_new <- B2b*b1_new
          colnames(B_new) <- rep(var_name2[1], 1)
          var_name3 <- var_name2[best.var]
          colnames(B_new) <- rep(var_name3, 1)
        }

        B_names <- paste(B3_names[best.var], B_name1, sep = "*")
        if (trunc.type == 2) B_names <- c(B_names, paste(B3_names[best.var], B_name2, sep = "*"))
        if (trunc.type == 1) B_names <- B_names

        var_name_list1 <- list()
        for (ll in 1:ncol(B_new)) {
          colnames(B_new)[ll] <- paste(var_name, colnames(B_new)[ll], sep = ":")
          var_name_list1 <- c(var_name_list1, list(colnames(B_new)[ll]))
          int.count1 <- int.count1 + 1
        }
      }

      if (int == FALSE) {
        var_name_list1 <- list()
        if (trunc.type == 2) {
          B_temp <- cbind(B, b1_new, b2_new) # Additive model with both truncated basis functions.
          B_new <- cbind(b1_new, b2_new)
          B_names <- c(B_name1, B_name2)
          var_name_list1 <- c(var_name_list1, list(var_name))
          var_name_list1 <- c(var_name_list1, list(var_name))  # Repeat it because there are two new truncated basis function in the set.
        }
        if (trunc.type == 1) {
          B_temp <- cbind(B, b1_new) # Additive model with one truncated basis function (i.e., the positive part).
          B_new <- b1_new
          B_names <- B_name1
          var_name_list1 <- c(var_name_list1, list(var_name))
        }
      }

      if (is.gee == TRUE | n > 1) meas_model <- score_fun_gee(Y, N, n_vec, VS.est_list, AWA.est_list, J2_list, Sigma2_list, J11.inv, JSigma11, mu.est, V.est, B_temp, B_new, nb, ...)
      if (is.gee == FALSE & n == 1) meas_model <- score_fun_glm(Y, N, VS.est_list, A_list, B1_list, mu.est, V.est, B_temp, B_new, nb, ...)

      score2 <- meas_model$score

      meas_model0 <- stat_out(Y, B_temp, TSS, GCV.null, pen, ...)
      GCVq2 <- meas_model0$GCVq1

      if (GCVq2 < (-10) | round(score2, 4) <= 0) {
        writeLines("** MARGE tolerance criteria met 1** \n")
        var.mod_temp <- c(var.mod_temp, NA)
        min_knot_vec_temp <- c(min_knot_vec_temp, NA)
        int.count1_temp <- c(int.count1_temp, NA)
        is.int_temp <- c(is.int_temp, int)

        trunc.type_temp <- c(trunc.type_temp, NA)
        X_red_temp <- c(X_red_temp, list(NA))
        B_new_list_temp <- c(B_new_list_temp, list(NA))
        var_name_list1_temp <- c(var_name_list1_temp, list(NA))
        B_names_temp <- c(B_names_temp, list(NA))
        B_temp_list <- c(B_temp_list, list(NA))

        score_term_temp <- c(score_term_temp, NA)

        if (length(var.mod_temp) == q) {
          breakFlag <- TRUE
          break
        }
        if (length(var.mod_temp) != q) next
      }

      if (GCVq2 >= (-10) | round(score2, 4) > 0) score_term_temp <- c(score_term_temp, score2)

      var.mod_temp <- c(var.mod_temp, var_name)
      min_knot_vec_temp <- c(min_knot_vec_temp, min_knot1)
      int.count1_temp <- c(int.count1_temp, int.count1)
      is.int_temp <- c(is.int_temp, int)
      trunc.type_temp <- c(trunc.type_temp, trunc.type)

      B_new_list_temp <- c(B_new_list_temp, list(B_new))
      var_name_list1_temp <- c(var_name_list1_temp, list(var_name_list1))
      B_names_temp <- c(B_names_temp, list(B_names))
      X_red_temp <- c(X_red_temp, list(X_red))
      B_temp_list <- c(B_temp_list, list(B_temp))

    }  # Terminate the for () loop to end v (variables) here.

    if (breakFlag == TRUE) break

    best.mod <- which.max(score_term_temp)  # Finds the best model (i.e., the max LOF) from candidate model/basis set. This becomes the new parent.

    score2 <- score_term_temp[best.mod]
    min_knot_vec1 <- min_knot_vec_temp[best.mod]
    int.count1 <- int.count1_temp[best.mod]
    int <- is.int_temp[best.mod]
    trunc.type <- trunc.type_temp[best.mod]

    B_new <- B_new_list_temp[[best.mod]]
    var_name_list1 <- var_name_list1_temp[[best.mod]]
    B_names <- B_names_temp[[best.mod]]
    X_red <- X_red_temp[[best.mod]]
    B_temp <- B_temp_list[[best.mod]]

    score_term <- c(score_term, score2)
    min_knot_vec <- c(min_knot_vec, min_knot_vec1)
    pred.name_vec <- c(pred.name_vec, colnames(B_new)[1])
    cut_vec <- c(cut_vec, X_red[min_knot_vec1])
    trunc.type_vec <- c(trunc.type_vec, trunc.type)
    is.int_vec <- c(is.int_vec, int)

    if (score_term[k + 1] < tols_score) {
      if (print.disp == TRUE) writeLines("\n ** MARGE tolerance criteria met 2** \n")
      breakFlag <- TRUE
      break
    }

    if (score_term[k + 1] >= tols_score) {
      if (int == TRUE) {
        mod_struct <- c(mod_struct, rep(c(rep(2, int.count1/trunc.type)), trunc.type))
        int.count <- int.count + 1
      }
      if (int == FALSE) mod_struct <- c(mod_struct, rep(1, trunc.type))

      B <- B_temp
      var_name_vec <- c(var_name_vec, colnames(B_new))
      var_name_list <- c(var_name_list, var_name_list1)
      B_names_vec <- c(B_names_vec, B_names)
      k <- k + 1
      m <- m + 2
    }

    if (nrow(B) <= (ncol(B) + 2)) { # To avoid the p > N, issue!
      if (print.disp == TRUE) writeLines("\n ** Parameter dimension exceeds N for MARGE ** \n")
      ok <- FALSE
    }

    if (m >= M) { # If model exceeds no. of set terms, terminate it.
      if (print.disp == TRUE) writeLines("\n ** Exceeded max no. of set terms for MARGE ** \n")
      ok <- FALSE
    }

    int
    B_names_vec
    mod_struct
  }

  colnames(B) <- B_names_vec
  B2 <- B

  # Algorithm 3 (backward pass) as in Friedman (1991) but for GLM/GEE use WIC.

  # WIC_vec_2 = WIC_vec_log = NA
  # ED: removing the use of both penalties
  WIC_vec = NA

  if (is.gee == FALSE) {
    # if (nb == FALSE) full.fit <- stats::glm(Y ~ B - 1, family = family)
    # if (nb == TRUE) full.fit <- gamlss::gamlss(Y ~ B - 1, family = "NBI", trace = FALSE)
    if (nb == TRUE) {
      full.fit <- gamlss::gamlss(Y ~ B - 1, family = "NBI", trace = FALSE)
    } else {
      stats::glm(Y ~ B - 1, family = family)
    }
  } else {
    full.fit <- geepack::geeglm(Y ~ B - 1, id = id, family = family, corstr = corstr)
    # ED: we can just pass the family argument and avoid the switch below
  }

  # if (is.gee == TRUE) {
  #   if (family == "gaussian") full.fit <- geepack::geeglm(Y ~ B - 1, id = id, corstr = corstr)
  #   if (family != "gaussian") full.fit <- geepack::geeglm(Y ~ B - 1, id = id, family = family, corstr = corstr)
  # }

  full.wic <- 0

  B_new <- B
  # cnames_2 <- list(colnames(B_new))
  # cnames_log <- list(colnames(B_new))
  # ED: removing the use of both penalties
  cnames <- list(colnames(B_new))

  # wic_mat_2 <- matrix(NA, ncol = ncol(B), nrow = ncol(B))
  # wic_mat_log <- matrix(NA, ncol = ncol(B), nrow = ncol(B))
  # ED: removing the use of both penalties
  wic_mat <- matrix(NA, ncol = ncol(B), nrow = ncol(B))

  # colnames(wic_mat_2) = colnames(wic_mat_log) = colnames(B)
  # ED: removing the use of both penalties
  colnames(wic_mat) = colnames(B)

  # wic_mat_2 <- cbind(wic_mat_2, rep(NA, ncol(B)))
  # wic_mat_log <- cbind(wic_mat_log, rep(NA, ncol(B)))
  # ED: removing the use of both penalties
  wic_mat <- cbind(wic_mat, rep(NA, ncol(B)))

  # colnames(wic_mat_2)[(ncol(B) + 1)] = colnames(wic_mat_log)[(ncol(B) + 1)] = "Forward pass model"
  # ED: removing the use of both penalties
  colnames(wic_mat)[(ncol(B) + 1)] = "Forward pass model"

  # wic_mat_2[1, (ncol(B) + 1)] = wic_mat_log[1, (ncol(B) + 1)] = full.wic
  # ED: removing the use of both penalties
  wic_mat[1, (ncol(B) + 1)] = full.wic

  # wic1_2 <- backward_sel_WIC(Y, N, n, B_new, id, family, corstr, nb, is.gee, ...)
  # wic1_log <- backward_sel_WIC(Y, N, n, B_new, id, family, corstr, nb, is.gee, ...)
  # ED: removing the use of both penalties
  wic1 <- backward_sel_WIC(Y, N, n, B_new, id, family, corstr, nb, is.gee, ...)

  # wic_mat_2[2, 2:(length(wic1_2) + 1)] <- wic1_2
  # wic_mat_log[2, 2:(length(wic1_log) + 1)] <- wic1_log
  # ED: removing the use of both penalties
  wic_mat[2, 2:(length(wic1) + 1)] <- wic1

  # WIC_2 <- sum(apply(wic_mat_2[1:2, ], 1, min, na.rm = TRUE)) + 2*ncol(B_new)
  # WIC_log <- sum(apply(wic_mat_log[1:2, ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new)
  # ED: removing the use of both penalties
  WIC <- switch(pen,
                two = sum(apply(wic_mat[1:2, ], 1, min, na.rm = TRUE)) + 2*ncol(B_new),
                log = sum(apply(wic_mat[1:2, ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new)
  )


  # WIC_vec_2 <- c(WIC_vec_2, WIC_2)
  # WIC_vec_log <- c(WIC_vec_log, WIC_log)
  # ED: removing the use of both penalties
  WIC_vec <- c(WIC_vec, WIC)

  # variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
  # variable.lowest_log <- as.numeric(which(wic1_log == min(wic1_log, na.rm = TRUE))[1])
  # ED: removing the use of both penalties
  variable.lowest <- as.numeric(which(wic1 == min(wic1, na.rm = TRUE))[1])

  # var.low.vec_2 <- c(colnames(B_new)[variable.lowest_2 + 1])
  # var.low.vec_log <- c(colnames(B_new)[variable.lowest_log + 1])
  # ED: removing the use of both penalties
  var.low.vec <- c(colnames(B_new)[variable.lowest + 1])

  # B_new_2 <- as.matrix(B_new[, -(variable.lowest_2 + 1)])
  # B_new_log <- as.matrix(B_new[, -(variable.lowest_log + 1)])
  # ED: removing the use of both penalties
  B_new_0 <- as.matrix(B_new[, -(variable.lowest + 1)])

  # cnames_2 <- c(cnames_2, list(colnames(B_new_2)))
  # cnames_log <- c(cnames_log, list(colnames(B_new_log)))
  # ED: removing the use of both penalties
  cnames <- c(cnames, list(colnames(B_new_0)))


  for (i in 2:(ncol(B) - 1)) {
    if (i != (ncol(B) - 1)) {
      # wic1_2 <- backward_sel_WIC(Y, N, n, B_new_2, id, family, corstr, nb, is.gee, ...)
      # wic1_log <- backward_sel_WIC(Y, N, n, B_new_log, id, family, corstr, nb, is.gee, ...)
      # ED: removing the use of both penalties
      wic1 <- backward_sel_WIC(Y, N, n, B_new_0, id, family, corstr, nb, is.gee, ...)

      # wic_mat_2[(i + 1), colnames(B_new_2)[-1]] <- wic1_2
      # wic_mat_log[(i + 1), colnames(B_new_log)[-1]] <- wic1_log
      # ED: removing the use of both penalties
      wic_mat[(i + 1), colnames(B_new_0)[-1]] <- wic1

      # WIC_2 <- sum(apply(wic_mat_2[1:(i + 1), ], 1, min, na.rm = TRUE)) + 2*ncol(B_new_2)
      # WIC_log <- sum(apply(wic_mat_log[1:(i + 1), ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new_log)
      # ED: removing the use of both penalties
      WIC <- switch(pen,
                    two = sum(apply(wic_mat[1:(i + 1), ], 1, min, na.rm = TRUE)) + 2*ncol(B_new_0),
                    log = sum(apply(wic_mat[1:(i + 1), ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new_0)
      )

      # WIC_vec_2 <- c(WIC_vec_2, WIC_2)
      # WIC_vec_log <- c(WIC_vec_log, WIC_log)
      # ED: removing the use of both penalties
      WIC_vec <- c(WIC_vec, WIC)

      # variable.lowest_2 <- as.numeric(which(wic1_2 == min(wic1_2, na.rm = TRUE))[1])
      # variable.lowest_log <- as.numeric(which(wic1_log == min(wic1_log, na.rm = TRUE))[1])
      # ED: removing the use of both penalties
      variable.lowest <- as.numeric(which(wic1 == min(wic1, na.rm = TRUE))[1])

      # var.low.vec_2 <- c(var.low.vec_2, colnames(B_new_2)[variable.lowest_2 + 1])
      # var.low.vec_log <- c(var.low.vec_log, colnames(B_new_log)[variable.lowest_log + 1])
      # ED: removing the use of both penalties
      var.low.vec <- c(var.low.vec, colnames(B_new_0)[variable.lowest + 1])

      # B_new_2 <- as.matrix(B_new_2[, -(variable.lowest_2 + 1)])
      # B_new_log <- as.matrix(B_new_log[, -(variable.lowest_log + 1)])
      # ED: removing the use of both penalties
      B_new_0 <- as.matrix(B_new_0[, -(variable.lowest + 1)])
    }

    if (i == (ncol(B) - 1)) {
      if (is.gee == FALSE) {
        if (nb == FALSE) {
          # full.fit_2 <- stats::glm(Y ~ B_new_2 - 1, family = family, ...)
          # full.fit_log <- stats::glm(Y ~ B_new_log - 1, family = family, ...)
          # ED: removing the use of both penalties
          full.fit <- stats::glm(Y ~ B_new_0 - 1, family = family, ...)

          # full.wald_2 <- (summary(full.fit_2)[12]$coef[-1, 3])^2
          # full.wald_log <- (summary(full.fit_log)[12]$coef[-1, 3])^2
          # ED: removing the use of both penalties
          full.wald <- (summary(full.fit)[12]$coef[-1, 3])^2

          # wic1_2 <- full.wald_2
          # wic1_log <- full.wald_log
          # ED: removing the use of both penalties
          wic1 <- full.wald
        }
        if (nb == TRUE) {
          # full.fit_2 <- gamlss::gamlss(Y ~ B_new_2 - 1, family = "NBI", trace = FALSE, ...)
          # full.fit_log <- gamlss::gamlss(Y ~ B_new_log - 1, family = "NBI", trace = FALSE, ...)
          # ED: removing the use of both penalties
          full.fit <- gamlss::gamlss(Y ~ B_new_0 - 1, family = "NBI", trace = FALSE, ...)

          sink(tempfile())
          # full.wald_2 <- ((as.matrix(summary(full.fit_2))[, 3])[-c(1, nrow(as.matrix(summary(full.fit_2))))])^2
          # full.wald_log <- ((as.matrix(summary(full.fit_log))[, 3])[-c(1, nrow(as.matrix(summary(full.fit_log))))])^2
          # ED: removing the use of both penalties
          full.wald <- ((as.matrix(summary(full.fit))[, 3])[-c(1, nrow(as.matrix(summary(full.fit))))])^2
          sink()

          # wic1_2 <- full.wald_2
          # wic1_log <- full.wald_log
          # ED: removing the use of both penalties
          wic1 <- full.wald
        }
      }

      if (is.gee == TRUE) {
        if (family == "gaussian") {
          # full.fit_2 <- geepack::geeglm(Y ~ B_new_2 - 1, id = id, corstr = corstr, ...)
          # full.fit_log <- geepack::geeglm(Y ~ B_new_log - 1, id = id, corstr = corstr, ...)
          # ED: removing the use of both penalties
          full.fit <- geepack::geeglm(Y ~ B_new_0 - 1, id = id, corstr = corstr, ...)

          # full.wald_2 <- (summary(full.fit_2)[6]$coef[-1, 3])^2
          # full.wald_log <- (summary(full.fit_log)[6]$coef[-1, 3])^2
          # ED: removing the use of both penalties
          full.wald <- (summary(full.fit)[6]$coef[-1, 3])^2
        }
        if (family != "gaussian") {
          # full.fit_2 <- geepack::geeglm(Y ~ B_new_2 - 1, id = id, family = family, corstr = corstr, ...)
          # full.fit_log <- geepack::geeglm(Y ~ B_new_log - 1, id = id, family = family, corstr = corstr, ...)
          # ED: removing the use of both penalties
          full.fit <- geepack::geeglm(Y ~ B_new_0 - 1, id = id, family = family, corstr = corstr, ...)

          # full.wald_2 <- (summary(full.fit_2)[6]$coef[-1, 3])^2
          # full.wald_log <- (summary(full.fit_log)[6]$coef[-1, 3])^2
          # ED: removing the use of both penalties
          full.wald <- (summary(full.fit)[6]$coef[-1, 3])^2
        }

        # wic1_2 <- full.wald_2
        # wic1_log <- full.wald_log
        # ED: removing the use of both penalties
        wic1 <- full.wald
      }

      # wic_mat_2[(i + 1), colnames(B_new_2)[-1]] <- wic1_2
      # wic_mat_log[(i + 1), colnames(B_new_log)[-1]] <- wic1_log
      # ED: removing the use of both penalties
      wic_mat[(i + 1), colnames(B_new_0)[-1]] <- wic1

      # WIC_2 <- sum(apply(wic_mat_2[1:(ncol(B)), ], 1, min, na.rm = TRUE)) + 2*ncol(B_new_2)
      # WIC_log <- sum(apply(wic_mat_log[1:(ncol(B)), ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new_log)
      # ED: removing the use of both penalties
      WIC <- switch(pen,
                    two = sum(apply(wic_mat[1:(ncol(B)), ], 1, min, na.rm = TRUE)) + 2*ncol(B_new_0),
                    log = sum(apply(wic_mat[1:(ncol(B)), ], 1, min, na.rm = TRUE)) + log(N)*ncol(B_new_0)
      )

      # WIC_vec_2 <- c(WIC_vec_2, WIC_2)
      # WIC_vec_log <- c(WIC_vec_log, WIC_log)
      # ED: removing the use of both penalties
      WIC_vec <- c(WIC_vec, WIC)

      # B_new_2 <- as.matrix(B_new_2[, -(variable.lowest_2)])
      # B_new_log <- as.matrix(B_new_log[, -(variable.lowest_log)])
      # ED: removing the use of both penalties
      B_new_0 <- as.matrix(B_new_0[, -(variable.lowest)])

      # colnames(B_new_2) <- "Intercept"
      # colnames(B_new_log) <- "Intercept"
      # ED: removing the use of both penalties
      colnames(B_new_0) <- "Intercept"
    }

    # cnames_2 <- c(cnames_2, list(colnames(B_new_2)))
    # cnames_log <- c(cnames_log, list(colnames(B_new_log)))
    # ED: removing the use of both penalties
    cnames <- c(cnames, list(colnames(B_new_0)))
  }

  # plot wic if required
  if (plot.wic) {
    # graphics::par(mfrow = c(1, 2), las = 0)
    # graphics::plot(rev(WIC_vec_2), type = "l", lwd = 2, ylab = "", xlab = "backward path", cex.lab = 2)
    # graphics::mtext(text = "WIC", line = 2, side = 2, cex = 2, las = 0)
    # graphics::title(main = bquote(paste( ~ paste(lambda) == .(2))), cex.main = 2.2, line = 1.4)
    # graphics::plot(rev(WIC_vec_log), type = "l", lwd = 2, ylab = "", xlab = "backward path", cex.lab = 2)
    # graphics::mtext(text = "WIC", line = 2, side = 2, cex = 2, las = 0)
    # graphics::title(main = bquote(paste( ~ paste(lambda) ~ " =  log(N)")), cex.main = 2.2, line = 1.4)
    # graphics::par(mfrow = c(1, 1))
    graphics::plot(rev(WIC_vec), type = "l", lwd = 2, ylab = "", xlab = "backward path", cex.lab = 2)
    graphics::mtext(text = "WIC", line = 2, side = 2, cex = 2, las = 0)
  }

  if (print.disp == TRUE) {
    writeLines("\n Forward pass output: \n")
    forw.info <- cbind(round(score_term, 4), pred.name_vec, cut_vec, trunc.type_vec, is.int_vec)
    colnames(forw.info) <- c("Score", "Predictor name", "Cut term (knot)", "No. of new parent terms", "Interaction?")
    print(forw.info)
  }

  # Some final model output, WIC, GCV etc.

  # ED: calculate the final B matrix (same for all model types)
  B_final <- as.matrix(B[, colnames(B)%in%cnames[[which.min(WIC_vec)]]])

  if (is.gee == TRUE) {
    # B_final <- as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]])
    # final_mod_2 <- geepack::geeglm(Y ~ B_final - 1, id = id, family = family, corstr = corstr, ...)
    #
    # B_final <- as.matrix(B[, colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]])
    # final_mod_log <- geepack::geeglm(Y ~ B_final - 1, id = id, family = family, corstr = corstr, ...)
    # ED: removing the use of both penalties
    final_mod <- geepack::geeglm(Y ~ B_final - 1, id = id, family = family, corstr = corstr, ...)

  }
  if (is.gee == FALSE) {
    if (nb == FALSE) {
      # B_final <- as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]])
      # final_mod_2 <- stats::glm(Y ~ B_final - 1, family = family, ...)
      # B_final <- as.matrix(B[, colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]])
      # final_mod_log <- stats::glm(Y ~ B_final - 1, family = family, ...)
      # ED: removing the use of both penalties
      final_mod <- stats::glm(Y ~ B_final - 1, family = family, ...)
    }
    if (nb == TRUE) {
      # B_final <- as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]])
      # final_mod.many_2 <- mvabund::manyglm(c(t(Y)) ~ B_final - 1, family = "negative.binomial", maxiter = 1000, maxiter2 = 100)
      # final_mod_2 <- MASS::glm.nb(c(t(Y)) ~ B_final - 1, method = "glm.fit2", init.theta = final_mod.many_2$theta, ...)
      #
      # B_final <- as.matrix(B[, colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]])
      # final_mod.many_log <- mvabund::manyglm(c(t(Y)) ~ B_final-1, family = "negative.binomial", maxiter = 1000, maxiter2 = 100)
      # final_mod_log <- MASS::glm.nb(c(t(Y)) ~ B_final - 1, method = "glm.fit2", init.theta = final_mod.many_log$theta, ...)
      # ED: removing the use of both penalties
      final_mod.many <- mvabund::manyglm(c(t(Y)) ~ B_final-1, family = "negative.binomial", maxiter = 1000, maxiter2 = 100)
      final_mod <- MASS::glm.nb(c(t(Y)) ~ B_final - 1, method = "glm.fit2", init.theta = final_mod.many$theta, ...)
    }
  }

  NN <- length(Y)

  # p_2 <- ncol(as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]))
  # df1a_2 <- p_2 + 2*(p_2 - 1)/2  # This matches the earth() package, SAS and Friedman (1991) penalty.
  # p_log <- ncol(as.matrix(B[, colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]))
  # df1a_log <- p_log + log(p_log - 1)/2  # This matches the earth() package, SAS and Friedman (1991) penalty.
  # ED: removing the use of both penalties
  p_0 <- ncol(as.matrix(B[, colnames(B)%in%cnames[[which.min(WIC_vec)]]]))
  df1a <- switch(pen,
                 two = p_0 + 2*(p_0 - 1)/2, # This matches the earth() package, SAS and Friedman (1991) penalty.
                 log = p_0 + log(p_log - 1)/2 # THIS NEEDS CHECKING
  )

  # RSS1_2 <- sum((Y - stats::fitted(final_mod_2))^2)
  # RSS1_log <- sum((Y - stats::fitted(final_mod_log))^2)
  # RSSq1_2 <- 1 - RSS1_2/sum((Y - mean(Y))^2)
  # RSSq1_log <- 1 - RSS1_log/sum((Y - mean(Y))^2)
  # GCV1_2 <- RSS1_2/(NN*(1 - df1a_2/NN)^2)
  # GCV1_log <- RSS1_log/(NN*(1 - df1a_log/NN)^2)
  # ED: removing the use of both penalties
  RSS1 <- sum((Y - stats::fitted(final_mod))^2)
  RSSq1 <- 1 - RSS1/sum((Y - mean(Y))^2)
  GCV1 <- RSS1/(NN*(1 - df1a/NN)^2)

  if (print.disp == TRUE) {
    # B_final <- as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]])
    # final_mod <- final_mod_2
    # wic_mat <- wic_mat_2
    # RSSq1 <- RSSq1_2
    # GCV1 <- GCV1_2
    # ED: the above is no longer needed

    writeLines("\n -- Final model (after pruning/backward selection) for MARGE -- \n")

    if (ncol(B_final) > 1) final_mat <- t(t(colnames(as.matrix(B_final))))
    if (ncol(B_final) == 1) final_mat <- as.matrix("Intercept", 1, 1)
    colnames(final_mat) <- "Selected variables in the final model:"
    print(final_mat)
    writeLines("\n Final model WIC: \n")
    print(min(wic_mat, na.rm = TRUE))
    writeLines("\n Final model Rsq: \n")
    print(RSSq1)
    writeLines("\n Final model GCV: \n")
    print(GCV1)
    writeLines("\n Final model coefs: \n")
    print(matrix(stats::coef(final_mod), dimnames = list(final_mat)))
  }

  # B_final <- list(as.matrix(B[, colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]), as.matrix(B[, colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]))
  # ED: no longer needed

  # wic_mat <- list(wic_mat_2, wic_mat_log)
  # min_wic_own <- list(min(wic_mat_2, na.rm = TRUE), (min(wic_mat_log, na.rm = TRUE)))
  # min_wic_own <-wic_mat
  # GCV1 <- list(GCV1_2, GCV1_log)
  # y_pred <- list(stats::predict(final_mod_2), stats::predict(final_mod_log))
  # final_mod <- list(final_mod_2, final_mod_log)
  y_pred <- stats::predict(final_mod)

  z <- NULL
  z$bx <- B_final
  z$wic_mat <- wic_mat
  z$min_wic_own <- min(wic_mat, na.rm = TRUE)
  z$GCV <- GCV1
  z$y_pred <- y_pred
  z$final_mod <- final_mod
  ## ED: adding some info required for predict.marge (and in turn, plotmo)
  z$coefficients <- final_mod$coefficients
  z$data <- data
  z$terms <- terms(formula)
  z$call <- match.call()
  z$y <- Y
  z$is.gee <- is.gee
  z$pen <- pen

  class(z) <- "marge"

  return(z)
}

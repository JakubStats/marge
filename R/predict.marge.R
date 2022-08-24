#' predict.marge
#'
#' @description S3 generic for predict.marge function for fitted MARS and MARGE objects when using the \code{marge} package.
#'
#' @param object : the final selected MARS/MARGE model. Only works for \code{mars_ls} and \code{marge} model objects.
#' @param newdata : the new set of predictor variables to make predictions on (test data). Must contain all predictors supplied within the original model formula.
#' @param is.marge : a logical argument, is this a MARGE model? The default is \code{FALSE}.
#' @param pen : the penalty used for the MARGE model (only applicable if \code{marge} was used). The default is \code{pen = "2"}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{predict.marge} returns a list of calculated values consisting of:
#' @return \code{eta.p} : the fitted linear predictor using the new (test) data.
#' @return \code{basis_new} : the model matrix for the new (test) data.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, 245--253.
#' @exportS3Method stats::predict marge
#' @importFrom stringr str_extract_all
#' @importFrom stats binomial poisson
#' @examples # Load the "leptrine" presence-absence data.
#'
#' data(leptrine)
#'
#' dat1 <- leptrine[[1]]   # Training data.
#' dat1_t <- leptrine[[2]] # Test data.
#'
#' Y <- dat1$Y             # Response variable.
#' N <- length(Y)          # Sample size (number of clusters).
#' n <- 1                  # Cluster size.
#' id <- rep(1:N, each = n)   # The ID of each cluster.
#'
#' X_pred <- dat1[, -c(3:10)]    # Design matrix using only two (of nine) predictors.
#' X_predt <- dat1_t[, -c(4:11)]
#'
#' # Set MARGE tuning parameters.
#'
#' family <- "binomial"    # The selected "exponential" family for the GLM/GEE.
#' is.gee <- FALSE         # Is the model a GEE?
#' nb <- FALSE             # Is this a negative binomial model?
#' tols_score <- 0.0001    # A set tolerance (stopping condition) in forward pass for MARGE.
#' M <- 21                 # A set threshold for the maximum number of basis functions to be used.
#' print.disp <- FALSE     # Print ALL the output?
#' pen <- 2                # Penalty to be used in GCV.
#' minspan <- NULL         # A set minimum span value.
#'
#' # Fit the MARGE models (about ~ 30 secs.)
#'
#' model_marge <- marge(X_pred, Y, N, n, id, family, corstr, pen, tols_score,
#'                      M, minspan, print.disp, nb, is.gee)
#'
#' # Predict on training data.
#'
#' pred_marge_2_y <- predict(model_marge, X_predt, X_pred, TRUE, "2")
#'
#' pred_marge_log_y <- predict(model_marge, X_predt, X_pred, TRUE, "logN")
# predict.marge <- function(object, newdata, X_pred, is.marge = FALSE, pen = "2", ...) {
#   if (is.marge == TRUE) {
#     if (pen == "2") mod <- object$final_mod[[1]]
#     if (pen == "logN") mod <- object$final_mod[[2]]
#   }
#
#   if (is.marge == FALSE) mod <- object$final_mod
#
#   if (ncol(stats::model.matrix(mod)) == 1) basis_new <- stats::model.matrix(mod)
#
#   if (ncol(stats::model.matrix(mod)) != 1) {
#     newdat_var <- colnames(newdata)[-1]  # Remove intercept term (1st column).
#
#     fitted_dat <- round(X_pred, digits = 4) # Original input data (need the variable names).
#
#     # Extract the variable names and cuts from the final model basis.
#
#     b1 <- substring(colnames(stats::model.matrix(mod)), 8)[-1]
#
#     var_name_list <- strsplit(b1, split = "\\*")
#     qq1 <- length(var_name_list)
#
#     pp1 <- length(newdat_var)
#     pp2 <- length(b1)
#
#     basis_new <- rep(1, nrow(newdata))
#
#     for (ww in 1:qq1) {
#       var1 <- var_name_list[[ww]]
#       qq2 <- length(var1)
#
#       if (qq2 == 1) {  # For additive structures.
#         var.vec_new1 <- c()
#
#         for (vv in 1:pp1) {
#           var.vec_new1 <- c(var.vec_new1, unique(unlist(stringr::str_extract_all(var1, newdat_var[vv]))))
#         }
#
#         var_name <- unlist(stringr::str_extract_all(var1, max(var.vec_new1)))
#
#         lenv1 <- length(strsplit(var1, "-|\\s")[[1]])
#
#         cut1 <- c(getNumberPart(var1))
#
#         if (lenv1 == 2) cut1 <- abs(cut1)
#         if (lenv1 == 3) cut1 <- cut1
#
#         if (length(cut1) > 1 & cut1[1] != cut1[2]) {
#           if (cut1[1] < 0 | cut1[2] < 0) {
#             var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#             cut00 <- which(cut1 < 0)
#             cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))
#
#             cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
#
#             if (floor(abs(cut1[cut00])) == floor(abs(cut3))) cut2 <- abs(cut1[cut00])*sign(cut3)
#           }
#           if (cut1[1] >= 0 & cut1[2] >= 0) {
#             var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#             cut00 <- which(cut1 == var_num0)
#             cut2 = cut3 = as.numeric(cut1[-cut00])
#           }
#         }
#         if (length(cut1) > 1 & (abs(cut1[1]) == abs(cut1[2]))) {
#           var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#           cut00 <- min(cut1)[1]
#           cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))
#
#           cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
#           cut2 <- abs(cut1[cut00])*sign(cut3)
#         }
#
#         if (length(cut1) == 1) cut2 = cut3 = cut1
#
#         temp_name <- paste("(", var_name, "-", cut2, ")", sep = "")
#         var_num <- which(colnames(newdata) == (noquote(var_name)))
#         var_chosen <- newdata[, var_num]
#
#         if (temp_name == var1) b1_new <- matrix(tp1(var_chosen, as.numeric(cut3)), ncol = 1)
#         if (temp_name != var1) b1_new <- matrix(tp2(var_chosen, as.numeric(cut3)), ncol = 1)
#
#         basis_new <- cbind(basis_new, round(b1_new, digits = 4))
#       }
#
#       if (qq2 == 2) {  # For interaction structures.
#         basis_new1 <- c()
#
#         for (yy in 1:qq2) {
#           var2 <- var1[yy]
#
#           var.vec_new1 <- c()
#
#           for (vv in 1:pp1) {
#             var.vec_new1 <- c(var.vec_new1, unique(unlist(stringr::str_extract_all(var2, newdat_var[vv]))))
#           }
#
#           var_name <- unlist(stringr::str_extract_all(var2, max(var.vec_new1)))
#
#           lenv1 <- length(strsplit(var2, "-|\\s")[[1]])
#
#           cut1 <- c(getNumberPart(var2))
#
#           if (lenv1 == 2) cut1 <- abs(cut1)
#           if (lenv1 == 3) cut1 <- cut1
#
#           if (length(cut1) > 1 & (cut1[1] != cut1[2])) {
#             if (cut1[1] < 0 | cut1[2] < 0) {
#               var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#               cut00 <- which(cut1 < 0)
#               cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))
#
#               cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
#
#               if (floor(abs(cut1[cut00])) == floor(abs(cut3))) cut2 <- abs(cut1[cut00])*sign(cut3)
#             }
#             if (cut1[1] >= 0 & cut1[2] >= 0) {
#               var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#               cut00 <- which(cut1 == var_num0)
#               cut2 = cut3 = as.numeric(cut1[-cut00])
#             }
#           }
#
#           if (length(cut1) > 1 & (abs(cut1[1]) == abs(cut1[2]))) {
#             var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
#             cut00 <- min(cut1)[1]
#             cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))
#
#             cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
#             cut2 <- abs(cut1[cut00])*sign(cut3)
#           }
#
#           if (length(cut1) == 1) cut2 = cut3 = cut1
#
#           temp_name <- paste("(", var_name, "-", cut2, ")", sep = "")
#           var_num <- which(colnames(newdata) == (noquote(var_name)))
#           var_chosen <- newdata[, var_num]
#
#           if (temp_name == var2) b1_new <- matrix(tp1(var_chosen, as.numeric(cut3)), ncol = 1)
#           if (temp_name != var2) b1_new <- matrix(tp2(var_chosen, as.numeric(cut3)), ncol = 1)
#
#           basis_new1 <- cbind(basis_new1, round(b1_new, digits = 4))
#         }
#
#         basis_new <- cbind(basis_new, basis_new1[, 1]*basis_new1[, 2])
#       }
#     }
#   }
#
#   eta.p <- c(stats::coef(mod)%*%t(basis_new))
#
#   list(eta.p = eta.p, basis_new = basis_new)
# }
predict.marge <- function(object, newdata, is.marge = TRUE, pen = c("2", "logN"), ...) {
  pen <- match.arg(pen)
  # obtain the appropriate model based on penalty (pen)
  if (is.marge) {
    mod <- switch(pen,
                  "2" = object$final_mod[[1]],
                  "logN" = object$final_mod[[2]])
  } else {
    mod <- object$final_mod
  }

  # assign "newdata" as the training data - will just produce fitted values if newdata is missing
  if (missing(newdata)) {
    newdata <- object$data
  }
  # obtain the variables
  newdat_var <- attr(object$terms, "term.labels")
  if (!all(newdat_var %in% colnames(newdata))) {
    stop("Original predictor terms from the model are not all found within 'newdata'")
  }

  if (length(mod$coefficients) == 1) {
    basis_new <- stats::model.matrix(mod)
  } else {

    # newdat_var <- colnames(newdata)[-1]  # Remove intercept term (1st column).

    fitted_dat <- round(object$data, digits = 4) # Original input data (need the variable names). ## ED why do we round here?

    # Extract the variable names and cuts from the final model basis.

    b1 <- substring(names(mod$coefficients), 8)[-1]

    var_name_list <- strsplit(b1, split = "\\*")
    qq1 <- length(var_name_list)

    pp1 <- length(newdat_var)
    pp2 <- length(b1)

    basis_new <- rep(1, nrow(newdata))

    for (ww in 1:qq1) {
      var1 <- var_name_list[[ww]]
      qq2 <- length(var1)

      if (qq2 == 1) {  # For additive structures.
        var.vec_new1 <- c()

        for (vv in 1:pp1) {
          var.vec_new1 <- c(var.vec_new1, unique(unlist(stringr::str_extract_all(var1, newdat_var[vv]))))
        }

        var_name <- unlist(stringr::str_extract_all(var1, max(var.vec_new1)))

        lenv1 <- length(strsplit(var1, "-|\\s")[[1]])

        cut1 <- c(getNumberPart(var1))

        if (lenv1 == 2) cut1 <- abs(cut1)
        if (lenv1 == 3) cut1 <- cut1

        if (length(cut1) > 1 & cut1[1] != cut1[2]) {
          if (cut1[1] < 0 | cut1[2] < 0) {
            var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
            cut00 <- which(cut1 < 0)
            cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))

            cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)

            if (floor(abs(cut1[cut00])) == floor(abs(cut3))) cut2 <- abs(cut1[cut00])*sign(cut3)
          }
          if (cut1[1] >= 0 & cut1[2] >= 0) {
            var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
            cut00 <- which(cut1 == var_num0)
            cut2 = cut3 = as.numeric(cut1[-cut00])
          }
        }
        if (length(cut1) > 1 & (abs(cut1[1]) == abs(cut1[2]))) {
          var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
          cut00 <- min(cut1)[1]
          cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))

          cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
          cut2 <- abs(cut1[cut00])*sign(cut3)
        }

        if (length(cut1) == 1) cut2 = cut3 = cut1

        temp_name <- paste("(", var_name, "-", cut2, ")", sep = "")
        var_num <- which(colnames(newdata) == (noquote(var_name)))
        var_chosen <- newdata[, var_num]

        if (temp_name == var1) b1_new <- matrix(tp1(var_chosen, as.numeric(cut3)), ncol = 1)
        if (temp_name != var1) b1_new <- matrix(tp2(var_chosen, as.numeric(cut3)), ncol = 1)

        basis_new <- cbind(basis_new, round(b1_new, digits = 4))
      }

      if (qq2 == 2) {  # For interaction structures.
        basis_new1 <- c()

        for (yy in 1:qq2) {
          var2 <- var1[yy]

          var.vec_new1 <- c()

          for (vv in 1:pp1) {
            var.vec_new1 <- c(var.vec_new1, unique(unlist(stringr::str_extract_all(var2, newdat_var[vv]))))
          }

          var_name <- unlist(stringr::str_extract_all(var2, max(var.vec_new1)))

          lenv1 <- length(strsplit(var2, "-|\\s")[[1]])

          cut1 <- c(getNumberPart(var2))

          if (lenv1 == 2) cut1 <- abs(cut1)
          if (lenv1 == 3) cut1 <- cut1

          if (length(cut1) > 1 & (cut1[1] != cut1[2])) {
            if (cut1[1] < 0 | cut1[2] < 0) {
              var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
              cut00 <- which(cut1 < 0)
              cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))

              cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)

              if (floor(abs(cut1[cut00])) == floor(abs(cut3))) cut2 <- abs(cut1[cut00])*sign(cut3)
            }
            if (cut1[1] >= 0 & cut1[2] >= 0) {
              var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
              cut00 <- which(cut1 == var_num0)
              cut2 = cut3 = as.numeric(cut1[-cut00])
            }
          }

          if (length(cut1) > 1 & (abs(cut1[1]) == abs(cut1[2]))) {
            var_num0 <- which(colnames(fitted_dat) == (noquote(var_name)))
            cut00 <- min(cut1)[1]
            cut01 <- abs(signif (as.numeric(cut1[cut00]), 4))

            cut3 <- utils::head(fitted_dat[which(abs(signif (fitted_dat[, var_num0], 4)) == cut01), var_num0], n = 1)
            cut2 <- abs(cut1[cut00])*sign(cut3)
          }

          if (length(cut1) == 1) cut2 = cut3 = cut1

          temp_name <- paste("(", var_name, "-", cut2, ")", sep = "")
          var_num <- which(colnames(newdata) == (noquote(var_name)))
          var_chosen <- newdata[, var_num]

          if (temp_name == var2) b1_new <- matrix(tp1(var_chosen, as.numeric(cut3)), ncol = 1)
          if (temp_name != var2) b1_new <- matrix(tp2(var_chosen, as.numeric(cut3)), ncol = 1)

          basis_new1 <- cbind(basis_new1, round(b1_new, digits = 4))
        }

        basis_new <- cbind(basis_new, basis_new1[, 1]*basis_new1[, 2])
      }
    }
  }

  eta.p <- c(stats::coef(mod)%*%t(basis_new))
  attr(eta.p, "basis_new") <- basis_new

  return(eta.p)
}

########################################################################
## MARGE_package.r
##
## Source file for the "marge" R-package.
########################################################################

## Required R-packages.

library(MASS);
library(geepack);
library(mvabund);
library(stringr);
library(gsubfn);
library(gamlss);
library(BayesTree);

#' @title Blue Mountains presence-absence data for the plant species \emph{Leptospermum trinervium}.
#' @description These data are binary (presence-absence) data collected on plant species for 8,678 sites located in the Blue Mountains region. It also contains environmental predictor variables collected at each site. There are are 1,751 absences and 6,927 presences. The `leptrine` data consists of training and test sets, such that \eqn{N} = 4,339. The nine environmental predictor variables given here have all been standardized.
#' @format A list containg two \code{matrices}: the training data with 4339 observations and 9 columns, and the test data with 4339 observations and 10 columns (it includes a column of zeroes for the intercept term). The columns are defined as follows:
#' \describe{
#'  \item{\code{RAIN_DRY_QTR}}{Recorded rainfall for the driest quater for each site.}
#'  \item{\code{FC}}{Recorded number of fire counts for each site.}
#'  \item{\code{TMP_MIN}}{Recorded minimum temperatures for each site.}
#'  \item{\code{TMP_SEAS}}{Recorded seasonal temperature for each site.}
#'  \item{\code{TMP_MN_WARM_QTR}}{Recorded minimum temperature for the warmest quater for each site.}
#'  \item{\code{RAIN_WET_QTR}}{Recorded rainfall for the westest quater for each site.}
#'  \item{\code{TMP_MN_COLD_QTR}}{Recorded minimum temperature for the coldest quater for each site.}
#'  \item{\code{TMP_MAX}}{Recorded maximum temperature for each site.}
#'  \item{\code{TMP_MN}}{Recorded minimum temperature for each site.}
#'  \item{\code{Y}}{Presence-absence for the plant species \emph{Leptospermum trinervium} for each site.}
#' }
#' @source The Blue Mountains presence-absence data for the speices \emph{Leptospermum trinervium} were obtained from \url{http://www.bionet.nsw.gov.au/}. Environmental data for Blue Mountains region: DRYAD entry doi:10.5061/dryad.985s5.
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @examples ## Load the data.
#'
#' data(leptrine)
#'
#' dat1<-leptrine[[1]]  # Training data.
#' Y<-dat1$Y            # Response variable.
#' N<-length(Y)         # Sample size (number of clusters).
#' n<-1                 # Cluster size.
#' id<-rep(1:N,each=n)  # The ID of each cluster.
#'
#' X_pred<-dat1[,-c(3:10)]  # Design matrix using only two (of nine) predictors.
#'
#' ## Set MARGE tunning parameters.
#'
#' family<-"binomial"   # The selected "exponential" family for the GLM/GEE.
#' is.gee<-FALSE        # Is the model a GEE?
#' nb<-FALSE            # Is this a negative binomial model?
#' tols_score<-0.0001   # A set tolerence (stopping condition) in forward pass for MARGE.
#' M<-21                # A set threshold for the maximum number of basis functions to be used.
#' print.disp<-FALSE    # Print ALL the output?
#' pen<-2               # Penalty to be used in GCV.
#' minspan<-NULL        # A set minimum span value.
#'
#' ## Fit the MARGE models (about ~ 30 secs.)
#'
#' mod<-marge(X_pred,Y,N,n,id,family,corstr,pen,tols_score,M,minspan,print.disp,nb,is.gee)
"leptrine"

#' tp1
#'
#' Truncated p-th power function (positive part).
#' @name tp1
#' @param x : a vector of predictor variable values.
#' @param t : a specified knot value.
#' @param p : the pth degree of the polynomial considered.
#' @return \code{tp1} returns a vector of values that have been transformed using a truncated p-th power function (positive part) for a specified knot value.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Stoklosa, J. and Warton, D.I. (2017). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, in review.
#' @export
#' @seealso \code{\link{tp2}}
#' @examples data(leptrine)
#'
#' dat1<-leptrine[[1]]
#' X_pred<-dat1[,1]   # One predicotr used.
#'
#' tp1(X_pred,1)      # Knot value set at x=1.
tp1<-function(x,t,p=1){((x-t)^p)*(x>t);}

#' tp2
#'
#' Truncated p-th power function (negative part).
#' @name tp2
#' @param x : a predictor variable value.
#' @param t : a specified knot value.
#' @param p : the pth degree of the polynomial considered.
#' @return \code{tp2} returns a vector of values that have been transformed using a truncated p-th power function (negative part) for a specified knot value.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{tp1}}
#' @examples data(leptrine)
#'
#' dat1<-leptrine[[1]]
#' X_pred<-dat1[,1]    # One predicotr used.
#'
#' tp2(X_pred,1)       # Knot value set at x=1.
tp2<-function(x,t,p=1){((t-x)^p)*(x<t);}

#' stat_out
#'
#' Fits a linear regresion model and calculates RSS/GCV measures (used for MARS linear models).
#' @name stat_out
#' @param Y : the response variable.
#' @param B1 : the model matrix of predictor variables.
#' @param TSS : total sum of squares.
#' @param GCV.null : GCV value for the intercept model.
#' @param pen : the set/fixed penalty used for the GCV (the default is 2).
#' @param ... : further arguments passed to or from other methods.
#' @details See the \code{earth} package for more details on the output measures calculated here.
#' @return \code{stat_out} returns a list of values, consisting of: RSS, RSSq1, GCV1 and GCVq1 values for the fitted model.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package=earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{stat_out_score_null}} and \code{\link{stat_out_score_glm_null}}
stat_out<-function(Y,B1,TSS,GCV.null,pen=2,...)
  {
  N<-length(Y);
  reg<-stats::lm.fit(B1,Y);

  if(sum(is.na(reg$coef))>0){RSS1=RSSq1=GCV1=GCVq1=NA;}
  if(sum(is.na(reg$coef))==0)
    {
    df1a<-ncol(B1)+pen*(ncol(B1)-1)/2;  # This matches the earth() package, SAS and Friedman (1991) penalty.

    RSS1<-sum((Y-stats::fitted(reg))^2);
    RSSq1<-1-RSS1/TSS;

    GCV1<-RSS1/(N*(1-(df1a)/N)^2);
    GCVq1<-1-GCV1/GCV.null;
    }

  list(RSS1=RSS1,RSSq1=round(RSSq1,10),GCV1=GCV1,GCVq1=round(GCVq1,10));
  }

#' stat_out_score_null
#'
#' A function that calculates parts for the score statistic for GEEs (it is used for the full path for forward selection).
#' @name stat_out_score_null
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param n : the cluster size.
#' @param id : the ID for each individual in the cluster.
#' @param family : the specified "exponential" family for GLMs. The default is \code{family="gaussian"}.
#' @param corstr : the specified "working correlation" structure. The default is \code{corstr="independence"}.
#' @param B_null : model matrix under the null model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param is.gee : a logical argument, is this a GEE model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @details The null model used here is by definition the current parent model. We compare the alternative model (a new basis function combined with the parent) with the null model. In the code used, only the null model is fit (the dispersion parameter is also estimated for GEE) and then the score statistic is obtained.
#' @return \code{stat_out_score_null} returns a list of values (mainly products of matrices) that make up the final score statistic calculation (required for another function).
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{stat_out}} and \code{\link{stat_out_score_glm_null}}
stat_out_score_null<-function(Y,N,n,id,family="gaussian",corstr="independence",B_null,nb=F,is.gee=F,...)
  {
  n_vec<-rep(n,N);

  if(is.gee==T)
    {
    if(family=="gaussian")
      {
      ests<-geepack::geeglm(Y~B_null-1,id=id,corstr=corstr);
      alpha.est<-ests$geese$alpha;
      }
    if(family!="gaussian")
      {
      ests<-geepack::geeglm(Y~B_null-1,id=id,family=family,corstr=corstr);
      alpha.est<-ests$geese$alpha;
      }

    mu.est<-as.matrix(stats::fitted.values(ests));
    }

  if(is.gee==F)
    {
    if(nb==T)
      {
      ests<-gamlss::gamlss(Y~B_null-1,family="NBI",trace=FALSE);
      mu.est<-as.matrix(stats::fitted.values(ests));
      }
    if(nb==F)
      {
      if(family=="gaussian"){ests<-stats::glm.fit(B_null,Y);}
      if(family=="binomial"){ests<-stats::glm.fit(B_null,Y,family=binomial(link = "logit"));}
      if(family=="poisson"){ests<-stats::glm.fit(B_null,Y,family=poisson(link = "log"));}
      mu.est<-as.matrix(stats::fitted.values(ests));
      }
    alpha.est<-1;
    }

  if(family=="gaussian"){V.est<-rep(1,nrow(mu.est));}
  if(family=="binomial"){V.est<-mu.est*(1-mu.est);}
  if(family=="poisson")
    {
    if(nb==F){V.est<-mu.est;}
    if(nb==T){V.est<-mu.est*(1+mu.est*(exp(ests$sigma.coef)));}
    }

  p<-ncol(B_null);
  n_vec1<-c(0,n_vec);

  VS.est_list<-list();
  AWA.est_list<-list();
  J2_list<-list();
  Sigma2_list<-list();

  J11<-matrix(0,nrow=p,ncol=p);
  Sigma11<-matrix(0,nrow=p,ncol=p);

  for(i in 1:N)
    {
    k<-sum(n_vec[1:i]);
    if(corstr=="independence"){R_alpha<-diag(1,nrow=n_vec[i],ncol=n_vec[i]);}
    if(corstr=="exchangeable"){R_alpha<-matrix(c(rep(alpha.est,n_vec[i]*n_vec[i])),ncol=n_vec[i])+diag(c(1-alpha.est),ncol=n_vec[i],nrow=n_vec[i]);}
    if(corstr=="ar1"){R_alpha<-alpha.est^outer(1:n_vec[i],1:n_vec[i],function(x,y){abs(x-y);});}

    V.est_i<-diag(sqrt(V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%R_alpha%*%diag(sqrt(V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i]);
    V.est_i_inv<-chol2inv(chol(V.est_i));

    S.est_i<-c(t(Y))[(sum(n_vec1[1:i])+1):k]-mu.est[(sum(n_vec1[1:i])+1):k];

    AWA.est_i<-V.est_i_inv%*%(S.est_i%*%t(S.est_i))%*%V.est_i_inv;

    if(nb==F){D.est_i<-diag((V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%B_null[(sum(n_vec1[1:i])+1):k,];}
    if(nb==T){D.est_i<-diag((mu.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%B_null[(sum(n_vec1[1:i])+1):k,];}

    J1_i<-t(D.est_i)%*%V.est_i_inv%*%D.est_i;
    J11<-J11+J1_i;
    J2_i<-t(D.est_i)%*%V.est_i_inv;

    Sigma1_i<-t(D.est_i)%*%AWA.est_i%*%(D.est_i);
    Sigma11<-Sigma11+Sigma1_i;
    Sigma2_i<-t(D.est_i)%*%AWA.est_i;

    VS.est_list<-c(VS.est_list,list(V.est_i_inv%*%S.est_i));
    AWA.est_list<-c(AWA.est_list,list(AWA.est_i));

    J2_list<-c(J2_list,list(J2_i));
    Sigma2_list<-c(Sigma2_list,list(Sigma2_i));
    }

  J11.inv<-chol2inv(chol(J11));

  JSigma11<-J11.inv%*%Sigma11%*%J11.inv;

  list(VS.est_list=VS.est_list,AWA.est_list=AWA.est_list,J2_list=J2_list,J11.inv=J11.inv,Sigma2_list=Sigma2_list,JSigma11=JSigma11,mu.est=mu.est,V.est=V.est);
  }

#' stat_out_score_glm_null
#'
#' @description A function that calculates parts of the score statistic for GLMs only (it is used for the full path for forward selection).
#' @name stat_out_score_glm_null
#' @param Y : the response variable.
#' @param family : the specified "exponential" family for the GLM. The default is \code{family="gaussian"}.
#' @param B_null : model matrix under the null model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{stat_out_score_glm_null} returns a list of values (mainly products of matrices) that make up the final score statistic calculation (required for another function).
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{stat_out}} and \code{\link{stat_out_score_glm_null}}
stat_out_score_glm_null<-function(Y,family="gaussian",B_null,nb=F,...)
  {
  if(nb==T)
    {
    ests<-gamlss::gamlss(Y~B_null-1,family="NBI",trace=FALSE);
    mu.est<-as.matrix(stats::fitted.values(ests));
    }
  if(nb==F)
    {
    if(family=="gaussian"){ests<-stats::glm.fit(B_null,Y);}
    if(family=="binomial"){ests<-stats::glm.fit(B_null,Y,family=binomial(link="logit"));}
    if(family=="poisson"){ests<-stats::glm.fit(B_null,Y,family=poisson(link="log"));}
    mu.est<-as.matrix(stats::fitted.values(ests));
    }

  if(family=="gaussian"){V.est<-rep(1,nrow(mu.est));}
  if(family=="binomial"){V.est<-mu.est*(1-mu.est);}
  if(family=="poisson")
    {
    if(nb==F){V.est<-mu.est;}
    if(nb==T){V.est<-mu.est*(1+mu.est*(exp(ests$sigma.coef)));}
    }

  VS.est_list<-(c(Y)-c(mu.est))/V.est;

  if(nb==F)
    {
    A_list<-chol2inv(chol((t(B_null)%*%diag(c(V.est))%*%B_null)));
    B1_list<-t(B_null)%*%diag(c(V.est));
    }

  if(nb==T)
    {
    A_list<-chol2inv(chol((t(B_null)%*%diag(c(mu.est^2/V.est))%*%B_null)));
    B1_list<-t(B_null)%*%diag(c(mu.est^2/V.est));
    }

  list(VS.est_list=VS.est_list,A_list=A_list,B1_list=B1_list,mu.est=mu.est,V.est=V.est);
  }

#' score_fun_glm
#'
#' Given estimates from the null model fit and the design matrix for alternative model, find the score statistic (this is used for GLMs only).
#' @name score_fun_glm
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param VS.est_list : a product of matrices.
#' @param A_list : a product of matrices.
#' @param B1_list : a product of matrices.
#' @param mu.est : esitmates of the fitted mean under the null model.
#' @param V.est : esitmates of the fitted variance under the null model.
#' @param B1 : model matrix under the null model.
#' @param XA : model matrix under the alternative model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{score_fun_glm} returns a calculated score statistic for the null and alternative model when fitting a GLM.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{score_fun_gee}}
score_fun_glm<-function(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B1,XA,nb=F,...)
  {
  reg<-try(stats::lm.fit(B1,Y),silent=TRUE);  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.

  if(class(reg)[1]=="try-error"){score<-NA;}

  if(class(reg)[1]!="try-error" && sum(is.na(reg$coef))>0){score<-NA;}

  if(class(reg)[1]!="try-error" && sum(is.na(reg$coef))==0)
    {
    VS.est_i<-unlist(VS.est_list);
    A_list_i<-unlist(A_list);
    B1_list_i<-unlist(B1_list);

    B_list_i<-B1_list_i%*%XA;

    if(nb==F)
      {
      D_list_i<-t(XA)%*%(XA*c(V.est));
      inv.XVX_22<-(D_list_i-t(B_list_i)%*%A_list_i%*%B_list_i);
      B.est<-t(V.est*VS.est_i)%*%XA;
      }
    if(nb==T)
      {
      D_list_i<-t(XA)%*%(XA*c((mu.est^2/V.est)));
      inv.XVX_22<-(D_list_i-t(B_list_i)%*%A_list_i%*%B_list_i);
      B.est<-t(((mu.est))*VS.est_i)%*%XA;
      }

    if(N<1500){score<-(B.est)%*%MASS::ginv(inv.XVX_22)%*%t(B.est);}
    if(N>1500){score<-(B.est)%*%chol2inv(chol(inv.XVX_22))%*%t(B.est);}
    }

  list(score=score);
  }

#' score_fun_gee
#'
#' Given estimates from the null and the design matrix from alternative model, find the score statistic (this is used for GEEs only).
#' @name score_fun_gee
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param n_vec : a vector of the number of cluster sizes.
#' @param VS.est_list : a product of matrices.
#' @param AWA.est_list : a product of matrices.
#' @param J2_list : a product of matrices.
#' @param Sigma2_list : a product of matrices.
#' @param J11.inv : a product of matrices.
#' @param JSigma11 : a product of matrices.
#' @param mu.est : esitmates of the fitted mean under the null model.
#' @param V.est : esitmates of the fitted variance under the null model.
#' @param B1 : model matrix under the null model.
#' @param XA : model matrix under the alternative model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{score_fun_gee} returns a calculated score statistic for the null and alternative model when fitting a GEE.
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{score_fun_glm}}
score_fun_gee<-function(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B1,XA,nb=F,...)
  {
  reg<-try(stats::lm.fit(B1,Y),silent=TRUE);  # This is not the model fit!! It just checks whether any issues occur for a simple linear regression model.

  if(class(reg)[1]=="try-error"){score<-NA;}

  if(class(reg)[1]!="try-error" && sum(is.na(reg$coef))>0){score<-NA;}

  if(class(reg)[1]!="try-error" && sum(is.na(reg$coef))==0)
    {
    p<-ncol(XA);
    p1<-ncol(B1)-ncol(XA);

    n_vec1<-c(0,n_vec);
    n<-max(n_vec);

    B.est<-matrix(0,nrow=p,ncol=1);
    Sigma22<-matrix(0,nrow=p,ncol=p);

    J21<-matrix(0,nrow=p,ncol=p1);
    Sigma21<-matrix(0,nrow=p,ncol=p1);

    for(i in 1:N)
      {
      k<-sum(n_vec[1:i]);

      VS.est_i<-VS.est_list[[i]];
      AWA.est_i<-AWA.est_list[[i]];
      J2_i<-J2_list[[i]];
      Sigma2_i<-Sigma2_list[[i]];

      if(nb==F){D.est_i<-diag((V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%XA[(sum(n_vec1[1:i])+1):k,];}
      if(nb==T){D.est_i<-diag((mu.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%XA[(sum(n_vec1[1:i])+1):k,];}

      J21<-J21+t(D.est_i)%*%t(J2_i);
      Sigma21<-Sigma21+t(D.est_i)%*%t(Sigma2_i);

      B.est<-B.est+t(D.est_i)%*%VS.est_i;
      Sigma22<-Sigma22+t(D.est_i)%*%AWA.est_i%*%(D.est_i);
      }

    Sigma<-Sigma22-(J21%*%J11.inv)%*%t(Sigma21)-(Sigma21%*%J11.inv)%*%t(J21)+J21%*%JSigma11%*%t(J21);

    score<-t(B.est)%*%MASS::ginv(Sigma)%*%B.est;
    }

  list(score=score);
  }

#' sand_fun
#'
#' A function that calculates the sandwich (or Pan's) standard error using estimates from the fitted GEE under an independent working correlation.
#' @name sand_fun
#' @param Y : the response variable.
#' @param X : the model matrix.
#' @param N : the number of clusters.
#' @param n_vec : a vector consisting of the number of cluster sizes.
#' @param mu.est : esitmates of the fitted mean function under the null model.
#' @param V.est : esitmates of the fitted variance function under the null model.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param omega : a regulariztion constant that is applied if your data is of high dimension (n>N). The default is \code{omega=1}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{sand_fun} returns the sandwich (or Pan's) standard error using estimates from the fitted GEE under an independent working correlation.
#' @details This function was used for the Arthropod example.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Liang, K. Y. and Zeger, S. L. (1986). Longitudinal aata analysis using generalized linear models. \emph{Biometrika}, 73, 13-22.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
sand_fun<-function(Y,X,N,n_vec,mu.est,V.est,nb=T,omega=1,...)
  {
  n<-max(n_vec);
  n_vec1<-c(0,n_vec);
  p<-ncol(X);

  A.est<-matrix(0,nrow=p,ncol=p);
  C.est<-matrix(0,nrow=p,ncol=p);
  R_u<-matrix(0,nrow=n,ncol=n);

  for(i in 1:N)
    {
    k<-sum(n_vec[1:i]);
    V.est_i<-diag(sqrt(V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%diag(n)%*%diag(sqrt(V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i]);
    if(nb==FALSE){D.est_i<-diag((V.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%X[(sum(n_vec1[1:i])+1):k,];}
    if(nb==TRUE){D.est_i<-diag((mu.est[(sum(n_vec1[1:i])+1):k]),nrow=n_vec[i],ncol=n_vec[i])%*%X[(sum(n_vec1[1:i])+1):k,];}
    S.est_i<-Y[(sum(n_vec1[1:i])+1):k]-mu.est[(sum(n_vec1[1:i])+1):k];
    A.est1<-t(D.est_i)%*%chol2inv(chol(V.est_i))%*%D.est_i;
    C.est1<-t(D.est_i)%*%chol2inv(chol(V.est_i))%*%(S.est_i%*%t(S.est_i))%*%chol2inv(chol(V.est_i))%*%D.est_i;
    A.est<-A.est+A.est1;
    C.est<-C.est+C.est1;
    }

  ## Pan's Robust modification estimator.

  V_est<-MASS::ginv(A.est)%*%C.est%*%MASS::ginv(A.est);

  list(V_est=V_est);
  }

#' backward_sel
#'
#' Backward selection funciton for MARS that uses the generalized cross validation criterion (GCV).
#' @name backward_sel
#' @param Y : the response variable.
#' @param B_new : the model matrix.
#' @param pen : the set/fixed penalty used for the GCV. The default is 2.
#' @param GCV.null : GCV value for the intercept model. The default is 0.001.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{backward_sel} returns the GCV from the fitted model.
#' @author Jakub Stoklosa and David I. Warton
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package=earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{backward_sel_WIC}}
backward_sel<-function(Y,B_new,pen=2,GCV.null=0.001,...)
  {
  N<-length(Y);

  GCV1<-c();

  for(j in 1:(ncol(B_new)-1))
    {
    B_new1<-as.matrix(B_new[,-(j+1)]);
    mod1<-stats::lm.fit(B_new1,Y);
    RSS_back<-sum((Y-stats::fitted(mod1))^2);

    p<-ncol(B_new1)+pen*(ncol(B_new1)-1)/2;  # This matches the earth() package, SAS and Friedman (1991) penalty.

    GCV1<-c(GCV1,1-(RSS_back/(N*(1-p/N)^2))/GCV.null);
    }

  return(GCV1);
  }

#' backward_sel_WIC
#'
#' Backward selection funciton for MARGE - uses the Wald information criterion (WIC).
#' @name backward_sel_WIC
#' @param Y : the response variable.
#' @param N : the number of clusters.
#' @param n : the cluster size.
#' @param B_new : the model matrix.
#' @param id : the ID for each individual in the cluster.
#' @param family : the specified "exponential" family for the GLM. The default is \code{family="gaussian"}.
#' @param corstr :the specified "working correlation" structure. The default is \code{corstr="independence"}.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param is.gee : is this a GEE model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{backward_sel_WIC} returns the Wald statistic from the fitted model (the penalty is apllied later on).
#' @author Jakub Stoklosa and David I. Warton
#' @references Stoklosa, J. Gibb, H. Warton, D.I. Fast forward selection for Generalized Estimating Equations With a Large Number of Predictor Variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{backward_sel}}
backward_sel_WIC<-function(Y,N,n,B_new,id,family="gaussian",corstr="independence",nb=F,is.gee=F,...)
  {
  if(is.gee==F)
    {
    if(nb==F)
      {
      fit<-stats::glm(c(t(Y))~B_new-1,family=family);
      wald1<-(summary(fit)[12]$coef[-1,3])^2;
      }
    if(nb==T)
      {
      sink(tempfile());
      fit<-gamlss::gamlss(Y~B_new-1,family="NBI",trace=FALSE);
      if(n==1){wald1<-((as.matrix(summary(fit))[,3])[-c(1,nrow(as.matrix(summary(fit))))])^2;}
      if(n>1)  # Used for the Arthropod data example.
        {
        n_vec<-rep(n,N);
        mu.est<-as.matrix(stats::fitted.values(fit));
        V.est<-mu.est*(1+mu.est*(exp(fit$sigma.coef)));
        wald1<-(fit$mu.coefficients[-1]/c(sqrt(diag(sand_fun(Y,B_new,N,n_vec,mu.est,V.est,T)$V_est)))[-1])^2;
        }
      sink();
      }
    }

  if(is.gee==T)
    {
    if(family=="gaussian")
      {
      fit<-geepack::geeglm(Y~B_new-1,id=id,corstr=corstr);
      wald1<-(summary(fit)[6]$coef[-1,3])^2;
      }
    if(family!="gaussian")
      {
      fit<-geepack::geeglm(Y~B_new-1,id=id,family=family,corstr=corstr);
      wald1<-(summary(fit)[6]$coef[-1,3])^2;
      }
    }

  wald1;
  }

#' max_span
#'
#' Truncates the predictor variable value to exclude extreme values in knots selection.
#' @name max_span
#' @param X_pred : a vector of values for the predictor variable.
#' @param q : the number of predictors used.
#' @param alpha : see Friedman (1991) equation (45). The default is 0.05.
#' @details Note that this equation comes from Friedman (1991) equation (45).
#' @return \code{max_span} returns a vector of trauncated predictor variable values.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
max_span<-function(X_pred,q,alpha=0.05)
  {
  N<-length(unique(X_pred));
  x<-sort(unique(X_pred));

  maxspan<-round((3-log2(alpha/q)));

  x_new<-x[-c(1:maxspan,floor(N-maxspan+1):N)];

  if(length(x_new)==0){return(x);}
  if(length(x_new)>0){return(x_new);}
  }

#' min_span
#'
#' A truncation function applied on the predictor variable for knot selection.
#' @name min_span
#' @param X_red : a vector of reduced predictor variable values.
#' @param q : the number of predictor variables used.
#' @param minspan : the set minimum span value. The default is \code{minspan=NULL}.
#' @param alpha : see Friedman (1991) equation (43). The default is 0.05.
#' @details  This function selects a minimum span between the knots to mitigate runs of correlated noise in the input data and hence avoiding estimation issues, this equation comes from Friedman (1991) equation 43.
#' @return \code{min_span} returns a vector of trauncated predictor variable values.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
min_span<-function(X_red,q,minspan=NULL,alpha=0.05)
  {
  N<-length((X_red));
  x<-sort((X_red));

  if(is.null(minspan)){minspan<-round((-log2(-(1/(q*N))*log(1-alpha))/2.5));}
  if(!is.null(minspan)){minspan<-minspan;}

  okA<-T;
  x_new<-min(x,na.rm=T);
  cc<-1;

  while(okA)
    {
    if((cc+minspan)>length(x))
      {
      break;
      okA<-F;
      }
    x_new1<-x[cc+(minspan+1)];
    x_new<-c(x_new,x_new1);
    cc<-cc+(minspan+1);
    }

  unique(x_new);
  }

#' mars_ls
#'
#' MARS estimation function for independent Gaussian response data. If non-Gaussian response data is used and a family function is specified (i.e., a GLM is used) then \code{mars_ls} applies a \code{glm} on the final selected basis function.
#' @name mars_ls
#' @param X_pred : a matrix of the predictor variables (note that this matrix excludes the intercept term).
#' @param Y : the response variable.
#' @param pen : a set penalty that is used for the GCV. The default is 2.
#' @param tols : a set tolerance for monitoring the convergence for the RSS between the parent and candidate model (this is the lack-of-fit criterion used for MARS). The default is 0.00001.
#' @param M : a set threshold for the number of basis functions to be used. The default is 21.
#' @param minspan : a set minimum span value. The default is \code{minspan=NULL}.
#' @param print.disp : a logical argument, do you want to print the output? The default is \code{FALSE}.
#' @param family : the specified "exponential" family if a GLM is used (only applied on the final model fit). The default is \code{family="gaussian"}.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param ... : further arguments passed to or from other methods.
#' @details Note that no argument is provided for an "additive structure only" model but one could just use the \code{gam()} function from \code{mgcv} instead.
#' @return \code{mars_ls} returns a list of calculated values consisting of:
#' @return \code{B_final} : the final selected basis model matrix.
#' @return \code{GCV_mat} : a matrix of GCV values for each fitted model.
#' @return \code{min_GCV_own} : the GCV for the final selected model.
#' @return \code{y_pred} : the fitted values from the final selected model.
#' @return \code{final_mod} : the final selected model matrix.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Milborrow, S. (2017a). Notes on the \code{earth} package. Package vignette. Available at: \url{http://127.0.0.1:31355/library/earth/doc/earth-notes.pdf}.
#' @references Milborrow, S. (2017b). \code{earth}: Multivariate Adaptive Regression Splines. R package version 4.4.7. Available at \url{http://CRAN.R-project.org/package=earth.}
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{marge}}
#' @examples ## Load the "leptrine" presence-absence data.
#'
#' data(leptrine)
#'
#' dat1<-leptrine[[1]]  # Training data.
#' Y<-dat1$Y            # Response variable.
#'
#' X_pred<-dat1[,-c(3:10)]  # Design matrix using only two (of nine) predictors.
#'
#' ## Set MARS tunning parameters.
#'
#' family<-"binomial"   # The selected "exponential" family.
#' nb<-FALSE            # Is this a negative binomial model?
#' tols<-0.0001         # A set tolerence (stopping condition) in the forward pass for MARS.
#' M<-21                # A set threshold for the maximum number of basis functions to be used.
#' print.disp<-FALSE    # Print ALL the output?
#' pen<-2               # Penalty to be used in GCV.
#' minspan<-NULL        # A set minimum span value.
#'
#' ## Fit the MARS models (about ~ 20 secs.)
#'
#' mod<-mars_ls(X_pred,Y,pen,tols,M,minspan,print.disp,family,nb)
mars_ls<-function(X_pred,Y,pen=2,tols=0.00001,M=21,minspan=NULL,print.disp=F,family="gaussian",nb=F,...)
  {
  N<-length(Y);    # Sample size.
  q<-ncol(X_pred);  # No. of predictor variables.

  ## Algorithm 2 (forward pass) as in Friedman (1991). Build an overfitted model using Rsq as the lack-of-fit (LOF).

  B<-as.matrix(rep(1,N)); # Start with the intercept model.

  colnames(B)<-c("Intercept");
  var_name_vec<-c("Intercept");
  var_name_list<-list("Intercept");
  B_names_vec<-c("Intercept");
  min_knot_vec<-c("Intercept");
  pred.name_vec<-c("Intercept");
  cut_vec<-c("Intercept");
  trunc.type_vec<-c(1);
  is.int_vec<-c("FALSE");
  mod_struct<-c(1);          # Univariate (1) or interaction (2).

  ## Null model setup - find the RSS, GCV and df. So initialize everything!

  B_null<-as.matrix(rep(1,N));
  RSS_term<-c(sum((Y-stats::fitted(stats::lm.fit(B_null,Y)))^2));
  TSS<-sum((Y-mean(Y))^2);
  RSSq_term<-c(1-RSS_term[1]/TSS);
  df1<-1;                     # Set DF to 1.
  GCV_term<-c(RSS_term/(N*(1-(df1)/N)^2));
  GCVq_term<-c(0);
  GCV.null<-GCV_term[1];

  m<-1;
  k<-1;
  breakFlag<-F;
  ok<-T;
  int.count<-0;

  ## Starting with loop 1 of algorithm considers terms (M). If model exceeds no. of set M we then terminate the loop.

  while(ok)
    {
    if(breakFlag==T){break;}

    var.mod_temp<-c();
    RSSq_term_temp<-c();
    GCVq_term_temp<-c();
    min_knot_vec_temp<-c();
    int.count1_temp<-c();
    is.int_temp<-c();
    trunc.type_temp<-c();

    B_new_list_temp<-list();
    var_name_list1_temp<-list();
    B_names_temp<-list();
    X_red_temp<-list();
    B_temp_list<-list();

    ## Loop 2 of algorithm considers all predictor variables in the data (q). Start with first variable.
    ## Note that variables can be repeated in the model, but cannot interact with itself (it doesn't
    ## really makes sense for it too anyway)!

    for(v in 1:q)
      {
      var_name<-colnames(X_pred)[v];
      X<-round(X_pred[,v],4);           # Round to 4 digits (to match earth()). X starts the new basis.

      ## These next few bits truncate the knot selection space: (1) reduce the space between knots; and (2) chop off ends bits.

      X_red1<-min_span(c(round(X,4)),q,minspan); # Reduce the space between knots.
      X_red2<-max_span(c(round(X,4)),q); # Truncate the ends of data to avoid extreme values.

      X_red<-intersect(X_red1,X_red2);

      RSSq_knot_both_int_mat<-c(); RSSq_knot_both_add_mat<-c();
      RSSq_knot_one_int_mat<-c(); RSSq_knot_one_add_mat<-c();
      GCVq_knot_both_int_mat<-c(); GCVq_knot_both_add_mat<-c();
      GCVq_knot_one_int_mat<-c(); GCVq_knot_one_add_mat<-c();

      int.count1<-0;

      ## This next bit is an indicator to check if the new variable can interact with any other variable already in the set. 0=no and >0, possible.

      if(ncol(B)>1){in.set<-sum(!var_name_vec%in%var_name);}
      if(ncol(B)==1){in.set<-0;}

      ## Loop 3 of algorithm considers ALL the knot locations (after trim) of the chosen variable in the current loop. The stratergy is simple for
      ## each knot find the lack-of-fit measures, and compare these with seperate fitted additive and interactions models.

      for(t in 1:length(X_red))
        {
        b1_new<-matrix(tp1(X,X_red[t]),ncol=1);  # Pairs of truncated power basis functions, positive and negative functions (first or both must be considered).
        b2_new<-matrix(tp2(X,X_red[t]),ncol=1);

        RSSq_knot_both_int<-c(); RSSq_knot_both_add<-c();
        RSSq_knot_one_int<-c(); RSSq_knot_one_add<-c();
        GCVq_knot_both_int<-c(); GCVq_knot_both_add<-c();
        GCVq_knot_one_int<-c(); GCVq_knot_one_add<-c();

        ## Find statistics for additive models only here, in fact this is only used once after initial entry of first variable, and used again if first variable was chosen.

        if(in.set==0)
          {
          B_new_both_add<-cbind(B,b1_new,b2_new);   # Additive model with both truncated functions.
          B_new_one_add<-cbind(B,b1_new);         # Additive model with one truncated function (positive part).

          meas_model_both_add<-stat_out(Y,B_new_both_add,TSS,GCV.null,pen); # Calculates the required stats.
          meas_model_one_add<-stat_out(Y,B_new_one_add,TSS,GCV.null,pen);

          RSSq_knot_both_add<-c(RSSq_knot_both_add,meas_model_both_add$RSSq1);
          RSSq_knot_one_add<-c(RSSq_knot_one_add,meas_model_one_add$RSSq1);
          GCVq_knot_both_add<-c(GCVq_knot_both_add,meas_model_both_add$GCVq1);
          GCVq_knot_one_add<-c(GCVq_knot_one_add,meas_model_one_add$GCVq1);

          RSSq_knot_both_int=RSSq_knot_one_int<--10000;      # Since no other vars. in the set then no interaction possible. Let the measures be some huge negative number.
          GCVq_knot_both_int=GCVq_knot_one_int<--10000;
          }

        ## Find LOF measures for all basis with interactions. Need to also fit additive models (they are candidates too)!

        if(in.set>0)
          {
          var_name_struct<-which(((var_name!=var_name_vec)*mod_struct)==1);  # This ensures that an interaction with the same variable is not included in the new basis.
          colnames(B)[1]<-c("");
          B2<-as.matrix(B[,var_name_struct]);
          if(k!=1 && (sum(!var_name_vec[-1]%in%var_name)>0)){B2<-as.matrix(B2[,-1]);}
          if(ncol(B2)==0){B2<-as.matrix(B[,1]);}

          for(nn in 1:ncol(B2))
            {
            B2a<-matrix(rep(B2[,nn],2),ncol=2);  # This is a basis function for potential parent basis. Need both part (hinges) here.
            B2b<-matrix(B2[,nn],ncol=1);
            B_new_both_int<-cbind(B,B2a*cbind(b1_new,b2_new));
            B_new_one_int<-cbind(B,B2b*b1_new);     # Interaction model with one truncated function (the positive part).

            meas_model_both_int<-stat_out(Y,B_new_both_int,TSS,GCV.null,pen);
            meas_model_one_int<-stat_out(Y,B_new_one_int,TSS,GCV.null,pen);

            RSSq_knot_both_int<-c(RSSq_knot_both_int,meas_model_both_int$RSSq1);
            RSSq_knot_one_int<-c(RSSq_knot_one_int,meas_model_one_int$RSSq1);
            GCVq_knot_both_int<-c(GCVq_knot_both_int,meas_model_both_int$GCVq1);
            GCVq_knot_one_int<-c(GCVq_knot_one_int,meas_model_one_int$GCVq1);
            }

          B_new_both_add<-cbind(B,b1_new,b2_new);
          B_new_one_add<-cbind(B,b1_new);
          meas_model_both_add<-stat_out(Y,B_new_both_add,TSS,GCV.null,pen);
          meas_model_one_add<-stat_out(Y,B_new_one_add,TSS,GCV.null,pen);
          RSSq_knot_both_add<-c(RSSq_knot_both_add,meas_model_both_add$RSSq1);
          RSSq_knot_one_add<-c(RSSq_knot_one_add,meas_model_one_add$RSSq1);
          GCVq_knot_both_add<-c(GCVq_knot_both_add,meas_model_both_add$GCVq1);
          GCVq_knot_one_add<-c(GCVq_knot_one_add,meas_model_one_add$GCVq1);
          }

        ## This next bit combines all lack-of-fit (LOF) measures that were found for each knot (t).

        RSSq_knot_both_int_mat<-rbind(RSSq_knot_both_int_mat,RSSq_knot_both_int);
        RSSq_knot_both_add_mat<-rbind(RSSq_knot_both_add_mat,RSSq_knot_both_add);
        RSSq_knot_one_int_mat<-rbind(RSSq_knot_one_int_mat,RSSq_knot_one_int);
        RSSq_knot_one_add_mat<-rbind(RSSq_knot_one_add_mat,RSSq_knot_one_add);
        GCVq_knot_both_int_mat<-rbind(GCVq_knot_both_int_mat,GCVq_knot_both_int);
        GCVq_knot_both_add_mat<-rbind(GCVq_knot_both_add_mat,GCVq_knot_both_add);
        GCVq_knot_one_int_mat<-rbind(GCVq_knot_one_int_mat,GCVq_knot_one_int);
        GCVq_knot_one_add_mat<-rbind(GCVq_knot_one_add_mat,GCVq_knot_one_add);

        } # Loop for t ends here.

      ## These four conditions are used to determine whether we include interactions (first) or additive parts only to the model, Then they look further by checking
      ## which trunctions to include for each parent basis in the set - i.e., one truncated function or both. All must be tested/checked for FUBR (NA) results!!!
      ## We apply the parsimonious principle here in case of tied measures - i.e., choose: additive model > interaction model, and one variable model > two variable model.

      if(sum(!(apply(RSSq_knot_both_int_mat,1,is.na)))==0 && sum(!(apply(RSSq_knot_one_int_mat,1,is.na)))==0)  # Look at both truncation types (check for NA). This says interactions were FUBR because of an NA present in both.
        {
        int<-F;
        if(sum(!is.na(RSSq_knot_both_add_mat))>0 && sum(!is.na(RSSq_knot_one_add_mat))>0)   # This says that an additive model for both one and two trunc types were OK.
          {
          if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
            {
            trunc.type<-2;
            RSSq_knot<-RSSq_knot_both_add_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            }
          if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
            {
            trunc.type<-1;
            RSSq_knot<-RSSq_knot_one_add_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            }
          }
        if(sum(!is.na(RSSq_knot_both_add_mat))==0 && sum(!is.na(RSSq_knot_one_add_mat))>0)   # This says additive for one trunc type was FUBR.
          {
          trunc.type<-1;           # For example, additive and both gave FUBR resutls, so use only the one (+) truncated function in the set.
          RSSq_knot<-RSSq_knot_one_add_mat;
          min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);  # We look for max because recall: GCVq1<-1-GCV1/GCV.null and RSSq1<-1-RSS1/TSS. Want these to be max. Tail used in case of ties.
          }
        if(sum(!is.na(RSSq_knot_both_add_mat))>0 && sum(!is.na(RSSq_knot_one_add_mat))==0) # This says additive for both trunc type was OK.
          {
          trunc.type<-2;
          RSSq_knot<-RSSq_knot_one_add_mat;
          min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
          }
        if(sum(!is.na(RSSq_knot_both_add_mat))==0 && sum(!is.na(RSSq_knot_one_add_mat))==0)  # This says everything (both additive trunc type) were FUBR.
          {
          breakFlag<-T;
          break;
          }
        }

      if(sum(!(apply(RSSq_knot_both_int_mat,1,is.na)))==0 && sum(!(apply(RSSq_knot_one_int_mat,1,is.na)))>0)  # This says interactions when using both trunc type were FUBR.
        {
        if(sum(!is.na(RSSq_knot_both_add_mat))==0)
          {
          trunc.type<-1;
          if(sum(!is.na(RSSq_knot_one_add_mat))==0)
            {
            int<-T;
            RSSq_knot<-RSSq_knot_one_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1)) # This says interactions using one trunc type were superior over the additive model.
            {
            int<-T;
            RSSq_knot<-RSSq_knot_one_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
            {
            int<-F;
            RSSq_knot<-RSSq_knot_one_add_mat;
            min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
            }
          }
        if(sum(!is.na(RSSq_knot_both_add_mat))>0)
          {
          if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1))
            {
            trunc.type<-1;
            if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-T;
              RSSq_knot<-RSSq_knot_one_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T))<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T)))
              {
              int<-F;
              RSSq_knot<-RSSq_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
              }
            }
          if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1))
            {
            int<-F;
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-2;
              RSSq_knot<-RSSq_knot_both_add_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              }
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              RSSq_knot<-RSSq_knot_one_add_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              }
            }
          }
        }

      if(sum(!(apply(RSSq_knot_both_int_mat,1,is.na)))>0 && sum(!(apply(RSSq_knot_one_int_mat,1,is.na)))==0) # This says interactions for one trunc. type was FUBR.
        {
        if(sum(!is.na(RSSq_knot_both_add_mat))==0)  # Check if additive for only one trunc type was OK.
          {
          if(sum(!is.na(RSSq_knot_one_add_mat))==0)
            {
            int<-T;
            trunc.type<-2;
            RSSq_knot<-RSSq_knot_both_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))  # Compare interaction (with both) with additive (with one).
            {
            int<-T;
            trunc.type<-2;
            RSSq_knot<-RSSq_knot_both_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
            {
            int<-F;
            trunc.type<-1;
            RSSq_knot<-RSSq_knot_one_add_mat;
            min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
            }
          }
        if(sum(!is.na(RSSq_knot_both_add_mat))>0)
          {
          if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1))
            {
            trunc.type<-1;
            if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-T;
              RSSq_knot<-RSSq_knot_both_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T))<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T)))
              {
              int<-F;
              RSSq_knot<-RSSq_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
              }
            }
          if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1))
            {
            int<-F;
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-2;
              RSSq_knot<-RSSq_knot_both_add_mat;
              min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
              }
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              RSSq_knot<-RSSq_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
              }
            }
          }
        }

      if(sum(!(apply(RSSq_knot_both_int_mat,1,is.na)))>0 && sum(!(apply(RSSq_knot_one_int_mat,1,is.na)))>0) # This says interactions for one trunc. and both type was OK.
        {
        if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
          {
          if(sum(!is.na(RSSq_knot_both_add_mat))>0)
            {
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>=utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-F;
              if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-2;
                RSSq_knot<-RSSq_knot_both_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
                }
              if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-1;
                RSSq_knot<-RSSq_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)));
                }
              }
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1))
              {
              if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-T;
                trunc.type<-2;
                RSSq_knot<-RSSq_knot_both_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
                }
              if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-F;
                trunc.type<-1;
                RSSq_knot<-RSSq_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
                }
              }
            }
          if(sum(!is.na(RSSq_knot_both_add_mat))==0)
            {
            if(sum(!is.na(RSSq_knot_one_add_mat))==0)
              {
              int<-T;
              trunc.type<-2;
              RSSq_knot<-RSSq_knot_both_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1)>=utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-F;
              trunc.type<-1;
              RSSq_knot<-RSSq_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
              }
            if(utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1)<utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-T;
              trunc.type<-2;
              RSSq_knot<-RSSq_knot_both_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
              }
            }
          }
        if(utils::tail(max(RSSq_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
          {
          if(sum(!is.na(RSSq_knot_both_add_mat))>0)
            {
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>=utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
              {
              int<-F;
              if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-2;
                RSSq_knot<-RSSq_knot_both_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
                }
              if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-1;
                RSSq_knot<-RSSq_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
                }
              }
            if(utils::tail(max(RSSq_knot_both_add_mat,na.rm=T),n=1)<utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-T;
                RSSq_knot<-RSSq_knot_one_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
                }
              if(utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-F;
                RSSq_knot<-RSSq_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)),n=1);
                }
              }
            }
          if(sum(!is.na(RSSq_knot_both_add_mat))==0)
            {
            if(sum(!is.na(RSSq_knot_one_add_mat))==0)
              {
              int<-T;
              trunc.type<-1;
              RSSq_knot<-RSSq_knot_one_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
              }
            if(sum(!is.na(RSSq_knot_one_add_mat))>0)
              {
              trunc.type<-1;
              if(utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1)>=utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
                {
                int<-F;
                RSSq_knot<-RSSq_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(RSSq_knot,6)));
                }
              if(utils::tail(max(RSSq_knot_one_add_mat,na.rm=T),n=1)<utils::tail(max(RSSq_knot_one_int_mat,na.rm=T),n=1))
                {
                int<-T;
                RSSq_knot<-RSSq_knot_one_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(RSSq_knot,6),na.rm=T),n=1)==round(RSSq_knot,6),arr.ind=T),n=1)[2];
                }
              }
            }
          }
        }

      b1_new<-matrix(tp1(X,X_red[min_knot1]),ncol=1);
      b2_new<-matrix(tp2(X,X_red[min_knot1]),ncol=1);
      colnames(b1_new)<-var_name;
      colnames(b2_new)<-var_name;

      B_name1<-paste("(",var_name,"-",signif(X_red[min_knot1],4),")",sep="");
      B_name2<-paste("(",signif(X_red[min_knot1],4),"-",var_name,")",sep="");

      if(int==T)
        {
        mod_struct1<-which(mod_struct==1);
        colnames(B)[1]<-c("");
        var_name1<-which(var_name_vec!=var_name);
        if(int.count==0){var_name2<-var_name_vec[var_name1];}
        if(int.count>0){var_name2<-var_name_vec[var_name_struct];}
        var_name_struct<-mod_struct1[mod_struct1%in%var_name1];
        B2<-as.matrix(B[,var_name_struct]);
        B3_names<-B_names_vec[var_name_struct];
        B3_names<-B3_names[-1];
        B2<-as.matrix(B2[,-1]);
        var_name2<-var_name2[-1];

        if(trunc.type==2)
          {
          B2a<-matrix(rep(B2[,best.var],2),ncol=2);
          B_temp<-cbind(B,B2a*cbind(b1_new,b2_new));
          B_new<-B2a*cbind(b1_new,b2_new);
          var_name3<-var_name2[best.var];
          colnames(B_new)<-rep(var_name3,2);
          }

        if(trunc.type==1)
          {
          B2b<-matrix(B2[,best.var],ncol=1);
          B_temp<-cbind(B,B2b*b1_new);     # Interaction model with one truncated function (i.e., the positive part).
          B_new<-B2b*b1_new;
          colnames(B_new)<-rep(var_name2[1],1);
          var_name3<-var_name2[best.var];
          colnames(B_new)<-rep(var_name3,1);
          }

        B_names<-paste(B3_names[best.var],B_name1,sep="*");
        if(trunc.type==2){B_names<-c(B_names,paste(B3_names[best.var],B_name2,sep="*"));}
        if(trunc.type==1){B_names<-B_names;}

        var_name_list1<-list();
        for(ll in 1:ncol(B_new))
          {
          colnames(B_new)[ll]<-paste(var_name,colnames(B_new)[ll],sep=":");
          var_name_list1<-c(var_name_list1,list(colnames(B_new)[ll]));
          int.count1<-int.count1+1;
          }

        pp<-ncol(B_temp);
        if(trunc.type==2){colnames(B_temp)[((pp-1):pp)]<-var_name_list1[[1]];}
        if(trunc.type==1){colnames(B_temp)[pp]<-var_name_list1[[1]];}
        }

      if(int==F)
        {
        var_name_list1<-list();
        if(trunc.type==2)
          {
          B_temp<-cbind(B,b1_new,b2_new); # Additive model with both truncated functions.
          B_new<-cbind(b1_new,b2_new);
          B_names<-c(B_name1,B_name2);
          var_name_list1<-c(var_name_list1,list(var_name));
          var_name_list1<-c(var_name_list1,list(var_name));  # Repeat it because there are new basis function selected.
          }
        if(trunc.type==1)
          {
          B_temp<-cbind(B,b1_new); # Additive model with one truncated function (positive part).
          B_new<-b1_new;
          B_names<-B_name1;
          var_name_list1<-c(var_name_list1,list(var_name));
          }
        }

      meas_model<-stat_out(Y,B_temp,TSS,GCV.null,pen);
      GCVq2<-meas_model$GCVq1;
      RSSq2<-meas_model$RSSq1;

      if(GCVq2<(-10) || (round(RSSq2,4)>(1-tols)))
        {
        writeLines("MARS Tolerance criteria met 1. \n");
        var.mod_temp<-c(var.mod_temp,NA);
        min_knot_vec_temp<-c(min_knot_vec_temp,NA);
        int.count1_temp<-c(int.count1_temp,NA);
        is.int_temp<-c(is.int_temp,int);

        trunc.type_temp<-c(trunc.type_temp,NA);
        X_red_temp<-c(X_red_temp,list(NA));
        B_new_list_temp<-c(B_new_list_temp,list(NA));
        var_name_list1_temp<-c(var_name_list1_temp,list(NA));
        B_names_temp<-c(B_names_temp,list(NA));
        B_temp_list<-c(B_temp_list,list(NA));

        RSSq_term_temp<-c(RSSq_term_temp,NA);
        GCVq_term_temp<-c(GCVq_term_temp,NA);

        if(length(var.mod_temp)==q)
          {
          breakFlag<-TRUE;
          break;
          }
        if(length(var.mod_temp)!=q){next;}
        }

      if(GCVq2>=(-10) || (round(RSSq2,4)<=(1-tols)))
        {
        RSSq_term_temp<-c(RSSq_term_temp,RSSq2);
        GCVq_term_temp<-c(GCVq_term_temp,GCVq2);
        }

      var.mod_temp<-c(var.mod_temp,var_name);
      min_knot_vec_temp<-c(min_knot_vec_temp,min_knot1);
      int.count1_temp<-c(int.count1_temp,int.count1);
      is.int_temp<-c(is.int_temp,int);
      trunc.type_temp<-c(trunc.type_temp,trunc.type);

      B_new_list_temp<-c(B_new_list_temp,list(B_new));
      var_name_list1_temp<-c(var_name_list1_temp,list(var_name_list1));
      B_names_temp<-c(B_names_temp,list(B_names));
      X_red_temp<-c(X_red_temp,list(X_red));
      B_temp_list<-c(B_temp_list,list(B_temp));

      } # End the for() v loop here.

    if(breakFlag==TRUE){break;}

    best.mod<-which.max(RSSq_term_temp);  # Finds the best model (max LOF) from your candidate model/basis set. This will become the new parent.

    RSSq2<-RSSq_term_temp[best.mod];
    GCVq2<-GCVq_term_temp[best.mod];
    min_knot_vec1<-min_knot_vec_temp[best.mod];
    int.count1<-int.count1_temp[best.mod];
    int<-is.int_temp[best.mod];
    trunc.type<-trunc.type_temp[best.mod];

    B_new<-B_new_list_temp[[best.mod]];
    var_name_list1<-var_name_list1_temp[[best.mod]];
    B_names<-B_names_temp[[best.mod]];
    X_red<-X_red_temp[[best.mod]];
    B_temp<-B_temp_list[[best.mod]];

    RSSq_term<-c(RSSq_term,RSSq2);
    GCVq_term<-c(GCVq_term,GCVq2);
    min_knot_vec<-c(min_knot_vec,min_knot_vec1);
    pred.name_vec<-c(pred.name_vec,colnames(B_new)[1]);
    cut_vec<-c(cut_vec,X_red[min_knot_vec1]);
    trunc.type_vec<-c(trunc.type_vec,trunc.type);
    is.int_vec<-c(is.int_vec,int);

    if(abs(RSSq_term[k]-RSSq_term[k+1])<tols)
      {
      if(print.disp==T){writeLines("\n ** MARS tolerance criteria met 2** \n");}
      breakFlag<-TRUE;
      break;
      }

    if(abs(RSSq_term[k]-RSSq_term[k+1])>=tols)
      {
      if(int==T)
        {
        mod_struct<-c(mod_struct,rep(c(rep(2,int.count1/trunc.type)),trunc.type));
        int.count<-int.count+1;
        }
      if(int==F){mod_struct<-c(mod_struct,rep(1,trunc.type));}
      B<-B_temp;
      var_name_vec<-c(var_name_vec,colnames(B_new));
      var_name_list<-c(var_name_list,var_name_list1);
      B_names_vec<-c(B_names_vec,B_names);
      k<-k+1;
      m<-m+2;
      }

    if(nrow(B)<=(ncol(B)+2)) # To avoid p>N issues!
      {
      if(print.disp==T){writeLines("\n ** MARS parameter dimension exceeds N ** \n");}
      ok<-F;
      }

    if(m>=M)   # If model exceeds no. of set terms, terminate it.
      {
      if(print.disp==T){writeLines("\n ** MARS exceeded max no. of set terms ** \n");}
      ok<-F;
      }

    int; B_names_vec; mod_struct; stat_out(Y,B,TSS,GCV.null,pen);
    }              # The term (m) for() loop ends here.

  colnames(B)<-B_names_vec;

  B2<-B; # Set this for final model to be the same as model selected from the forward pass.

  ## Algorithm 3 (backward pass), as in Friedman (1991). This uses GCV as the selection criterion.

  p<-ncol(B)+pen*(ncol(B)-1)/2;  # This matches the earth() package, SAS and Friedman (1991) penalty.

  full_RSS<-sum((Y-stats::fitted(stats::lm.fit(B-1,Y)))^2);
  full_GCV<-full_RSS/(N*(1-p/N)^2);
  full_GCV<-1-(full_GCV/GCV.null);

  B_new<-B;
  GCV_mat<-matrix(NA,ncol=ncol(B),nrow=ncol(B));
  colnames(GCV_mat)<-colnames(B);
  GCV_mat<-cbind(GCV_mat,rep(NA,ncol(B)));
  colnames(GCV_mat)[(ncol(B)+1)]<-"Forward pass model";

  GCV_mat[1,(ncol(B)+1)]<-full_GCV;
  GCV1<-backward_sel(Y,B_new,pen,GCV.null);
  GCV_mat[2,2:(length(GCV1)+1)]<-GCV1;
  variable.lowest<-utils::tail(which(GCV1==max(GCV1,na.rm=T)),n=1);
  var.low.vec<-c(colnames(B_new)[variable.lowest+1]);
  B_new<-as.matrix(B_new[,-(variable.lowest+1)]);

  for(i in 2:(ncol(B)-1))
    {
    GCV1<-backward_sel(Y,B_new,pen,GCV.null);

    variable.lowest<-utils::tail(which(GCV1==max(GCV1,na.rm=T)),n=1);
    var.low.vec<-c(var.low.vec,colnames(B_new)[variable.lowest+1]);

    if(i!=(ncol(B)-1)){GCV_mat[(i+1),colnames(B_new)[-1]]<-GCV1;}
    if(i==(ncol(B)-1)){GCV_mat[(i+1),1]<-GCV1;}

    B_new<-as.matrix(B_new[,-(variable.lowest+1)]);
    }

  number.rm1<-which(max(GCV_mat,na.rm=T)==GCV_mat,arr.ind=TRUE)[1]-1;

  if(number.rm1==0)
    {
    if(print.disp==T){writeLines("\n ** Forward pass model was chosen after pruning/backward selection for MARS** \n");}
    B_final<-B;
    }

  if(number.rm1==(nrow(GCV_mat)-1))
    {
    if(print.disp==T){writeLines("\n ** Intercept model was chosen after pruning/backward selection for MARS** \n");}
    B_final<-as.matrix(B[,1]);
    }

  if(!(number.rm1==0) && !(number.rm1==(nrow(GCV_mat)-1)))
    {
    number.rm2<-c();
    var.low.vec_red<-var.low.vec[1:number.rm1];
    for(k in 1:length(var.low.vec_red))
      {
      number.rm2<-c(number.rm2,which(var.low.vec_red[k]==colnames(B)));
      }
    B_final<-B[,-c(number.rm2)];
    }

  if(ncol(B_final)==1){colnames(B_final)<-"Intercept";}

  if(print.disp==T)
    {
    writeLines("\n Forward pass (MARS) output: \n");
    forw.info<-cbind(round(GCVq_term,4),round(RSSq_term,4),pred.name_vec,cut_vec,trunc.type_vec,is.int_vec);
    colnames(forw.info)<-c("GCVq","RSSq","Predictor name","Cut term (knot)","No. of new parent terms","Interaction?");
    print(forw.info);
    }

  ## Some measures to be reported, GCV, RSq, etc. for the final model fit.

  if(nb==F){final_mod<-stats::glm(Y~B_final-1,family=family);}

  if(nb==T)
    {
    est.many<-mvabund::manyglm(Y~B_final-1,family="negative.binomial",maxiter=1000,maxiter2=100);
    final_mod<-MASS::glm.nb(c(t(Y))~B_final-1,method="glm.fit2",init.theta=est.many$theta);
    }

  final_stats<-stat_out(Y,B_final,TSS,GCV.null,pen);

  if(print.disp==T)
    {
    writeLines("\n -- Final model was chosen (after pruning/backward selection) for MARS -- \n");
    final_mat<-t(t(colnames(as.matrix(B_final))));
    colnames(final_mat)<-"Selected variables in the final model:";
    print(final_mat);
    writeLines("\n Final model GCV: \n");

    print(final_stats$GCV1);
    writeLines("\n Final model Rsq: \n");
    print(final_stats$RSSq1);
    writeLines("\n Final model coefs: \n");
    print(matrix(stats::coef(final_mod),dimnames=list(final_mat)));
    }

  z<-NULL;
  z$bx<-B_final;
  z$GCV_mat<-GCV_mat;
  z$min_GCV_own<-final_stats$GCV1;
  z$y_pred<-stats::predict(final_mod);
  z$final_mod<-final_mod;

  class(z)<-"marge";
  return(z);
  }

#' marge
#'
#' MARS fitting function for generalized linear models (GLM) and generalized estimating equations (GEE).
#' @name marge
#' @param X_pred : a matrix of the predictor variables. Note that this matrix should include a column of 1's for the intercept term.
#' @param Y : the response variable. A vector of length n by N.
#' @param N : the number of clusters.
#' @param n : the cluster size.
#' @param id : a vector which identifies the clusters. The length of \code{id} should be the same as the number of observations. Data are assumed to be sorted so that observations on a cluster are contiguous rows for all entities in the formula.
#' @param family : the specified family for the GLM/GEE. The default is \code{family="gaussian"}. The current available families are: "gaussian", "binomial", "poisson" and the "negative binomial" (to get that use family="poisson" and nb=T).
#' @param corstr : the specified "working correlation" structure for the GEE. The default is \code{corstr="independence"}.
#' @param pen : a set penalty used for the GCV (note: MARGE doesn't actually use this). The default is 2.
#' @param tols_score : the set tolerance for monitoring the convergence for the difference in score statistics between the parent and candidate model (this is the lack-of-fit criterion used for MARGE). The default is 0.00001
#' @param M : a set threshold for the number of basis functions to be used. The default is 21.
#' @param minspan : a set minimum span value. The default is \code{minspan=NULL}.
#' @param print.disp : a logical argument, do you want to print the output? The default is \code{FALSE}.
#' @param nb : a logical argument, is the model a negative binomial model? The default is \code{FALSE}.
#' @param is.gee : is this a GEE model? The default is \code{FALSE}.
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
#' @references Friedman, J. (1991). Multivariate adaptive regression splines. \emph{The Annals of Statistics}, 19, 1-67.
#' @references Stoklosa, J., Gibb, H. and Warton, D.I. (2014). Fast forward selection for generalized estimating equations with a large number of predictor variables. \emph{Biometrics}, 70, 110-120.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @seealso \code{\link{mars_ls}} and \code{\link{backward_sel_WIC}}
#' @importFrom gamlss gamlss
#' @importFrom geepack geeglm
#' @importFrom mvabund manyglm
#' @importFrom MASS glm.nb
#' @examples ## Load the "leptrine" presence-absence data.
#'
#' data(leptrine)
#'
#' dat1<-leptrine[[1]]  # Training data.
#' Y<-dat1$Y            # Response variable.
#' N<-length(Y)         # Sample size (number of clusters).
#' n<-1                 # Cluster size.
#' id<-rep(1:N,each=n)  # The ID of each cluster.
#'
#' X_pred<-dat1[,-c(3:10)] # Design matrix using only two (of nine) predictors.
#'
#' ## Set MARGE tunning parameters.
#'
#' family<-"binomial"   # The selected "exponential" family for the GLM/GEE.
#' is.gee<-FALSE        # Is the model a GEE?
#' nb<-FALSE            # Is this a negative binomial model?
#' tols_score<-0.0001   # A set tolerence (stopping condition) in forward pass for MARGE.
#' M<-21                # A set threshold for the maximum number of basis functions to be used.
#' print.disp<-FALSE    # Print ALL the output?
#' pen<-2               # Penalty to be used in GCV.
#' minspan<-NULL        # A set minimum span value.
#'
#' ## Fit the MARGE models (about ~ 30 secs.)
#'
#' mod<-marge(X_pred,Y,N,n,id,family,corstr,pen,tols_score,M,minspan,print.disp,nb,is.gee)
marge<-function(X_pred,Y,N,n=1,id=c(1:length(Y)),family="gaussian",corstr="independence",pen=2,tols_score=0.00001,M=21,minspan=NULL,print.disp=F,nb=F,is.gee=F,...)
  {
  NN<-length(Y);    # Total sample size = N*n.

  q<-ncol(X_pred);  # No. of predictor variables.

  n_vec<-rep(n,N);

  ## Algorithm 2 (forward pass) as in Friedman (1991). Uses score statistics instead of RSS, etc.

  B<-as.matrix(rep(1,NN));       # Start with the intercept model.

  colnames(B)<-c("Intercept");
  var_name_vec<-c("Intercept");
  var_name_list<-list("Intercept");
  B_names_vec<-c("Intercept");
  min_knot_vec<-c("Intercept");
  pred.name_vec<-c("Intercept");
  cut_vec<-c("Intercept");
  trunc.type_vec<-c(1);
  is.int_vec<-c("FALSE");
  mod_struct<-c(1);          # Univariate (1) or interaction (2).

  score_term<-c(0);

  TSS<-sum((Y-mean(Y))^2);
  GCV.null<-TSS/(NN*(1-1/NN)^2);

  ## Null model setup.

  m<-1;
  k<-1;
  breakFlag<-F;
  ok<-T;
  int.count<-0;

  while(ok)
    {
    if(breakFlag==T){break;}

    var.mod_temp<-c();
    score_term_temp<-c();
    min_knot_vec_temp<-c();
    int.count1_temp<-c();
    is.int_temp<-c();
    trunc.type_temp<-c();

    B_new_list_temp<-list();
    var_name_list1_temp<-list();
    B_names_temp<-list();
    X_red_temp<-list();
    B_temp_list<-list();

    ## Obtain/calculate the null stats here (speeds things up).

    if(is.gee==T | n>1)
      {
      B_null_stats<-stat_out_score_null(Y,N,n,id,family,corstr,B,nb,is.gee);

      VS.est_list<-B_null_stats$VS.est_list;
      AWA.est_list<-B_null_stats$AWA.est_list;
      J2_list<-B_null_stats$J2_list;
      J11.inv<-B_null_stats$J11.inv;
      Sigma2_list<-B_null_stats$Sigma2_list;
      JSigma11<-B_null_stats$JSigma11;
      mu.est<-B_null_stats$mu.est;
      V.est<-B_null_stats$V.est;
      }
    if(is.gee==F & n==1)
      {
      B_null_stats<-stat_out_score_glm_null(Y,family,B,nb);

      VS.est_list<-B_null_stats$VS.est_list;
      A_list<-B_null_stats$A_list;
      B1_list<-B_null_stats$B1_list;
      mu.est<-B_null_stats$mu.est;
      V.est<-B_null_stats$V.est;
      }

    for(v in 1:q)
      {
      var_name<-colnames(X_pred)[v];
      X<-round(X_pred[,v],4);

      X_red1<-min_span(c(round(X,4)),q,minspan);  # Reduce the space between knots.
      X_red2<-max_span(c(round(X,4)),q);  # Truncate the ends of data to avoid extreme values.

      X_red<-intersect(X_red1,X_red2);

      score_knot_both_int_mat<-c(); score_knot_both_add_mat<-c();
      score_knot_one_int_mat<-c(); score_knot_one_add_mat<-c();

      int.count1<-0;

      if(ncol(B)>1){in.set<-sum(!var_name_vec%in%var_name)}
      if(ncol(B)==1){in.set<-0;}

      for(t in 1:length(X_red))
        {
        b1_new<-matrix(tp1(X,X_red[t]),ncol=1);  # Pairs of truncated functions.
        b2_new<-matrix(tp2(X,X_red[t]),ncol=1);

        score_knot_both_int<-c(); score_knot_both_add<-c();
        score_knot_one_int<-c(); score_knot_one_add<-c();

        if(in.set==0)
          {
          B_new_both_add<-cbind(B,b1_new,b2_new);   # Additive model with both truncated functions.
          B_new_one_add<-cbind(B,b1_new);         # Additive model with one truncated function (positive part).

          if(is.gee==T | n>1){meas_model_both_add<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_both_add,cbind(b1_new,b2_new),nb);}
          if(is.gee==F & n==1){meas_model_both_add<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_both_add,cbind(b1_new,b2_new),nb);}
          if(is.gee==T | n>1){meas_model_one_add<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_one_add,b1_new,nb);}
          if(is.gee==F & n==1){meas_model_one_add<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_one_add,b1_new,nb);}

          score_knot_both_add<-c(score_knot_both_add,meas_model_both_add$score);
          score_knot_one_add<-c(score_knot_one_add,meas_model_one_add$score);

          score_knot_both_int=score_knot_one_int= -100000;   # Interaction set is impossible since there is nothing to interact with, so let the LOF measure be a huge negative number.
          }

        if(in.set>0)
          {
          var_name_struct<-which(((var_name!=var_name_vec)*mod_struct)==1);
          colnames(B)[1]<-c("");
          B2<-as.matrix(B[,var_name_struct]);
          if(k!=1 && (sum(!var_name_vec[-1]%in%var_name)>0)){B2<-as.matrix(B2[,-1]);}
          if(ncol(B2)==0){B2<-as.matrix(B[,1]);}

          for(nn in 1:ncol(B2))
            {
            B2a<-matrix(rep(B2[,nn],2),ncol=2);
            B2b<-matrix(B2[,nn],ncol=1);
            B_new_both_int<-cbind(B,B2a*cbind(b1_new,b2_new));
            B_new_one_int<-cbind(B,B2b*b1_new);     # Interaction model with one truncated function (i.e., the positive part).

            if(is.gee==T | n>1){meas_model_both_int<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_both_int,B2a*cbind(b1_new,b2_new),nb);}
            if(is.gee==F & n==1){meas_model_both_int<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_both_int,B2a*cbind(b1_new,b2_new),nb);}
            if(is.gee==T | n>1){meas_model_one_int<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_one_int,B2b*b1_new,nb);}
            if(is.gee==F & n==1){meas_model_one_int<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_one_int,B2b*b1_new,nb);}

            score_knot_both_int<-c(score_knot_both_int,meas_model_both_int$score);
            score_knot_one_int<-c(score_knot_one_int,meas_model_one_int$score);
            }

          B_new_both_add<-cbind(B,b1_new,b2_new);
          B_new_one_add<-cbind(B,b1_new);

          if(is.gee==T | n>1){meas_model_both_add<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_both_add,cbind(b1_new,b2_new),nb);}
          if(is.gee==F & n==1){meas_model_both_add<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_both_add,cbind(b1_new,b2_new),nb);}
          if(is.gee==T | n>1){meas_model_one_add<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_new_one_add,b1_new,nb);}
          if(is.gee==F & n==1){meas_model_one_add<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_new_one_add,b1_new,nb);}

          score_knot_both_add<-c(score_knot_both_add,meas_model_both_add$score);
          score_knot_one_add<-c(score_knot_one_add,meas_model_one_add$score);
          }

        score_knot_both_int_mat<-rbind(score_knot_both_int_mat,score_knot_both_int);
        score_knot_both_add_mat<-rbind(score_knot_both_add_mat,score_knot_both_add);
        score_knot_one_int_mat<-rbind(score_knot_one_int_mat,score_knot_one_int);
        score_knot_one_add_mat<-rbind(score_knot_one_add_mat,score_knot_one_add);
        }

      ## See the LM code above in regards to what the conditions below actually do.

      if(sum(!(apply(score_knot_both_int_mat,1,is.na)))==0 && sum(!(apply(score_knot_one_int_mat,1,is.na)))==0)
        {
        int<-F;
        if(sum(!is.na(score_knot_both_add_mat))>0 && sum(!is.na(score_knot_one_add_mat))>0)
          {
          if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
            {
            trunc.type<-2;
            score_knot<-score_knot_both_add_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
            }
          if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
            {
            trunc.type<-1;
            score_knot<-score_knot_one_add_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
            }
          }
        if(sum(!is.na(score_knot_both_add_mat))==0 && sum(!is.na(score_knot_one_add_mat))>0)
          {
          trunc.type<-1;
          score_knot<-score_knot_one_add_mat;
          min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
          }
        if(sum(!is.na(score_knot_both_add_mat))>0 && sum(!is.na(score_knot_one_add_mat))==0)
          {
          trunc.type<-2;
          score_knot<-score_knot_one_add_mat;
          min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
          }
        if(sum(!is.na(score_knot_both_add_mat))==0 && sum(!is.na(score_knot_one_add_mat))==0)
          {
          breakFlag<-T;
          break;
          }
        }

      if(sum(!(apply(score_knot_both_int_mat,1,is.na)))==0 && sum(!(apply(score_knot_one_int_mat,1,is.na)))>0)
        {
        if(sum(!is.na(score_knot_both_add_mat))==0)
          {
          trunc.type<-1;
          if(sum(!is.na(score_knot_one_add_mat))==0)
            {
            int<-T;
            score_knot<-score_knot_one_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
            }
          if(sum(!is.na(score_knot_one_add_mat))>0)
            {
            if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-T;
              score_knot<-score_knot_one_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-F;
              score_knot<-score_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            }
          }
        if(sum(!is.na(score_knot_both_add_mat))>0)
          {
          if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1))
            {
            trunc.type<-1;
            if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-T;
              score_knot<-score_knot_one_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              int<-F;
              score_knot<-score_knot_one_add_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            }
          if(utils::tail(max(score_knot_one_int_mat,na.rm=T))<=utils::tail(max(score_knot_both_add_mat,na.rm=T)))
            {
            int<-F;
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              score_knot<-score_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-2;
              score_knot<-score_knot_both_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            }
          }
        }

      if(sum(!(apply(score_knot_both_int_mat,1,is.na)))>0 && sum(!(apply(score_knot_one_int_mat,1,is.na)))==0)
        {
        if(sum(!is.na(score_knot_both_add_mat))==0)
          {
          if(sum(!is.na(score_knot_one_add_mat))==0)
            {
            int<-T;
            trunc.type<-2;
            score_knot<-score_knot_both_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
            {
            int<-T;
            trunc.type<-2;
            score_knot<-score_knot_both_int_mat;
            min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
            best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
            }
          if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
            {
            int<-F;
            trunc.type<-1;
            score_knot<-score_knot_one_add_mat;
            min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
            }
          }
        if(sum(!is.na(score_knot_both_add_mat))>0)
          {
          if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1))
            {
            trunc.type<-2;
            if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1))
              {
              int<-T;
              score_knot<-score_knot_both_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1))
              {
              int<-F;
              score_knot<-score_knot_both_add_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              }
            }
          if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1))
            {
            int<-F;
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              score_knot<-score_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
              {
              trunc.type<-2;
              score_knot<-score_knot_both_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            }
          }
        }

      if(sum(!(apply(score_knot_both_int_mat,1,is.na)))>0 && sum(!(apply(score_knot_one_int_mat,1,is.na)))>0)
        {
        if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
          {
          if(sum(!is.na(score_knot_both_add_mat))>0)
            {
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>=utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-F;
              if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-2;
                score_knot<-score_knot_both_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
                }
              if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-1;
                score_knot<-score_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)));
                }
              }
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1))
              {
              if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-T;
                trunc.type<-2;
                score_knot<-score_knot_both_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
                }
              if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-F;
                trunc.type<-1;
                score_knot<-score_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
                }
              }
            }
          if(sum(!is.na(score_knot_both_add_mat))==0)
            {
            if(utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1)>=utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-F;
              trunc.type<-1;
              score_knot<-score_knot_one_add_mat;
              min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
              }
            if(utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1)<utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1))
              {
              int<-T;
              trunc.type<-2;
              score_knot<-score_knot_both_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            }
          }
        if(utils::tail(max(score_knot_both_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
          {
          if(sum(!is.na(score_knot_both_add_mat))>0)
            {
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>=utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
              {
              int<-F;
              if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-1;
                score_knot<-score_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
                }
              if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                trunc.type<-2;
                score_knot<-score_knot_both_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
                }
              }
            if(utils::tail(max(score_knot_both_add_mat,na.rm=T),n=1)<utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
              {
              trunc.type<-1;
              if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)>utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-T;
                score_knot<-score_knot_one_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
                }
              if(utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1)<=utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1))
                {
                int<-F;
                score_knot<-score_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)),n=1);
                }
              }
            }
          if(sum(!is.na(score_knot_both_add_mat))==0)
            {
            trunc.type<-1;
            if(sum(!is.na(score_knot_one_add_mat))==0)
              {
              int<-T;
              score_knot<-score_knot_one_int_mat;
              min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
              best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
              }
            if(sum(!is.na(score_knot_one_add_mat))>0)
              {
              if(utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1)>=utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
                {
                int<-F;
                score_knot<-score_knot_one_add_mat;
                min_knot1<-utils::tail(which.max(round(score_knot,6)));
                }
              if(utils::tail(max(score_knot_one_add_mat,na.rm=T),n=1)<utils::tail(max(score_knot_one_int_mat,na.rm=T),n=1))
                {
                int<-T;
                score_knot<-score_knot_one_int_mat;
                min_knot1<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[1];
                best.var<-utils::tail(which(utils::tail(max(round(score_knot,6),na.rm=T),n=1)==round(score_knot,6),arr.ind=T),n=1)[2];
                }
              }
            }
          }
        }

      b1_new<-matrix(tp1(X,X_red[min_knot1]),ncol=1);
      b2_new<-matrix(tp2(X,X_red[min_knot1]),ncol=1);
      colnames(b1_new)<-var_name;
      colnames(b2_new)<-var_name;

      B_name1<-paste("(",var_name,"-",signif(X_red[min_knot1],4),")",sep="");
      B_name2<-paste("(",signif(X_red[min_knot1],4),"-",var_name,")",sep="");

      if(int==T)
        {
        mod_struct1<-which(mod_struct==1);
        colnames(B)[1]<-c("");
        var_name1<-which(var_name_vec!=var_name);
        if(int.count==0){var_name2<-var_name_vec[var_name1];}
        if(int.count>0){var_name2<-var_name_vec[var_name_struct];}
        var_name_struct<-mod_struct1[mod_struct1%in%var_name1];
        B2<-as.matrix(B[,var_name_struct]);
        B3_names<-B_names_vec[var_name_struct];
        B3_names<-B3_names[-1];
        B2<-as.matrix(B2[,-1]);
        var_name2<-var_name2[-1];

        if(trunc.type==2)
          {
          B2a<-matrix(rep(B2[,best.var],2),ncol=2);
          B_temp<-cbind(B,B2a*cbind(b1_new,b2_new));
          B_new<-B2a*cbind(b1_new,b2_new);
          var_name3<-var_name2[best.var];
          colnames(B_new)<-rep(var_name3,2);
          }

        if(trunc.type==1)
          {
          B2b<-matrix(B2[,best.var],ncol=1);
          B_temp<-cbind(B,B2b*b1_new);     # Interaction model with one truncated basis function (i.e., the positive part).
          B_new<-B2b*b1_new;
          colnames(B_new)<-rep(var_name2[1],1);
          var_name3<-var_name2[best.var];
          colnames(B_new)<-rep(var_name3,1);
          }

        B_names<-paste(B3_names[best.var],B_name1,sep="*");
        if(trunc.type==2){B_names<-c(B_names,paste(B3_names[best.var],B_name2,sep="*"));}
        if(trunc.type==1){B_names<-B_names;}

        var_name_list1<-list();
        for(ll in 1:ncol(B_new))
          {
          colnames(B_new)[ll]<-paste(var_name,colnames(B_new)[ll],sep=":");
          var_name_list1<-c(var_name_list1,list(colnames(B_new)[ll]));
          int.count1<-int.count1+1;
          }
        }

      if(int==F)
        {
        var_name_list1<-list();
        if(trunc.type==2)
          {
          B_temp<-cbind(B,b1_new,b2_new); # Additive model with both truncated basis functions.
          B_new<-cbind(b1_new,b2_new);
          B_names<-c(B_name1,B_name2);
          var_name_list1<-c(var_name_list1,list(var_name));
          var_name_list1<-c(var_name_list1,list(var_name));  # Repeat it because there are two new truncated basis function in the set.
          }
        if(trunc.type==1)
          {
          B_temp<-cbind(B,b1_new); # Additive model with one truncated basis function (i.e., the positive part).
          B_new<-b1_new;
          B_names<-B_name1;
          var_name_list1<-c(var_name_list1,list(var_name));
          }
        }

      if(is.gee==T | n>1){meas_model<-score_fun_gee(Y,N,n_vec,VS.est_list,AWA.est_list,J2_list,Sigma2_list,J11.inv,JSigma11,mu.est,V.est,B_temp,B_new,nb);}
      if(is.gee==F & n==1){meas_model<-score_fun_glm(Y,N,VS.est_list,A_list,B1_list,mu.est,V.est,B_temp,B_new,nb);}

      score2<-meas_model$score;

      meas_model0<-stat_out(Y,B_temp,TSS,GCV.null,pen);
      GCVq2<-meas_model0$GCVq1;

      if(GCVq2<(-10) || round(score2,4)<=0)
        {
        writeLines("** MARGE tolerance criteria met 1** \n");
        var.mod_temp<-c(var.mod_temp,NA);
        min_knot_vec_temp<-c(min_knot_vec_temp,NA);
        int.count1_temp<-c(int.count1_temp,NA);
        is.int_temp<-c(is.int_temp,int);

        trunc.type_temp<-c(trunc.type_temp,NA);
        X_red_temp<-c(X_red_temp,list(NA));
        B_new_list_temp<-c(B_new_list_temp,list(NA));
        var_name_list1_temp<-c(var_name_list1_temp,list(NA));
        B_names_temp<-c(B_names_temp,list(NA));
        B_temp_list<-c(B_temp_list,list(NA));

        score_term_temp<-c(score_term_temp,NA);

        if(length(var.mod_temp)==q)
          {
          breakFlag<-TRUE;
          break;
          }
        if(length(var.mod_temp)!=q){next;}
        }

      if(GCVq2>=(-10) || round(score2,4)>0)
        {
        score_term_temp<-c(score_term_temp,score2);
        }

      var.mod_temp<-c(var.mod_temp,var_name);
      min_knot_vec_temp<-c(min_knot_vec_temp,min_knot1);
      int.count1_temp<-c(int.count1_temp,int.count1);
      is.int_temp<-c(is.int_temp,int);
      trunc.type_temp<-c(trunc.type_temp,trunc.type);

      B_new_list_temp<-c(B_new_list_temp,list(B_new));
      var_name_list1_temp<-c(var_name_list1_temp,list(var_name_list1));
      B_names_temp<-c(B_names_temp,list(B_names));
      X_red_temp<-c(X_red_temp,list(X_red));
      B_temp_list<-c(B_temp_list,list(B_temp));

      }  # Terminate the for() loop to end v (variables) here.

    if(breakFlag==TRUE){break;}

    best.mod<-which.max(score_term_temp);  # Finds the best model (i.e., the max LOF) from your candidate model/basis set. This becomes the new parent.

    score2<-score_term_temp[best.mod];
    min_knot_vec1<-min_knot_vec_temp[best.mod];
    int.count1<-int.count1_temp[best.mod];
    int<-is.int_temp[best.mod];
    trunc.type<-trunc.type_temp[best.mod];

    B_new<-B_new_list_temp[[best.mod]];
    var_name_list1<-var_name_list1_temp[[best.mod]];
    B_names<-B_names_temp[[best.mod]];
    X_red<-X_red_temp[[best.mod]];
    B_temp<-B_temp_list[[best.mod]];

    score_term<-c(score_term,score2);
    min_knot_vec<-c(min_knot_vec,min_knot_vec1);
    pred.name_vec<-c(pred.name_vec,colnames(B_new)[1]);
    cut_vec<-c(cut_vec,X_red[min_knot_vec1]);
    trunc.type_vec<-c(trunc.type_vec,trunc.type);
    is.int_vec<-c(is.int_vec,int);

    if(score_term[k+1]<tols_score)
      {
      if(print.disp==T){writeLines("\n ** MARGE tolerance criteria met 2** \n");}
      breakFlag<-TRUE;
      break;
      }

    if(score_term[k+1]>=tols_score)
      {
      if(int==T)
        {
        mod_struct<-c(mod_struct,rep(c(rep(2,int.count1/trunc.type)),trunc.type));
        int.count<-int.count+1;
        }
      if(int==F){mod_struct<-c(mod_struct,rep(1,trunc.type));}

      B<-B_temp;
      var_name_vec<-c(var_name_vec,colnames(B_new));
      var_name_list<-c(var_name_list,var_name_list1);
      B_names_vec<-c(B_names_vec,B_names);
      k<-k+1;
      m<-m+2;
      }

    if(nrow(B)<=(ncol(B)+2)) # To aviod the p>N, issue!
      {
      if(print.disp==T){writeLines("\n ** Parameter dimension exceeds N for MARGE ** \n");}
      ok<-F;
      }

    if(m>=M) # If model exceeds no. of set terms, terminate it.
      {
      if(print.disp==T){writeLines("\n ** Exceeded max no. of set terms for MARGE ** \n");}
      ok<-F;
      }

    int; B_names_vec; mod_struct;
    }

  colnames(B)<-B_names_vec;
  B2<-B;

  ## Algorithm 3 (backward pass) as in Friedman (1991) but for GLM/GEE use WIC.

  WIC_vec_2=WIC_vec_log=NA;

  if(is.gee==F)
    {
    if(nb==F){full.fit<-stats::glm(Y~B-1,family=family);}
    if(nb==T){full.fit<-gamlss::gamlss(Y~B-1,family="NBI",trace=FALSE);}
    }

  if(is.gee==T)
    {
    if(family=="gaussian"){full.fit<-geepack::geeglm(Y~B-1,id=id,corstr=corstr);}
    if(family!="gaussian"){full.fit<-geepack::geeglm(Y~B-1,id=id,family=family,corstr=corstr);}
    }

  full.wic<-0;

  B_new<-B;
  cnames_2<-list(colnames(B_new));
  cnames_log<-list(colnames(B_new));

  wic_mat_2<-matrix(NA,ncol=ncol(B),nrow=ncol(B));
  wic_mat_log<-matrix(NA,ncol=ncol(B),nrow=ncol(B));
  colnames(wic_mat_2)=colnames(wic_mat_log)=colnames(B);
  wic_mat_2<-cbind(wic_mat_2,rep(NA,ncol(B)));
  wic_mat_log<-cbind(wic_mat_log,rep(NA,ncol(B)));
  colnames(wic_mat_2)[(ncol(B)+1)]=colnames(wic_mat_log)[(ncol(B)+1)]="Forward pass model";

  wic_mat_2[1,(ncol(B)+1)]=wic_mat_log[1,(ncol(B)+1)]=full.wic;

  wic1_2<-backward_sel_WIC(Y,N,n,B_new,id,family,corstr,nb,is.gee);
  wic1_log<-backward_sel_WIC(Y,N,n,B_new,id,family,corstr,nb,is.gee);

  wic_mat_2[2,2:(length(wic1_2)+1)]<-wic1_2;
  wic_mat_log[2,2:(length(wic1_log)+1)]<-wic1_log;

  WIC_2<-sum(apply(wic_mat_2[1:2,],1,min,na.rm=T))+2*ncol(B_new);
  WIC_log<-sum(apply(wic_mat_log[1:2,],1,min,na.rm=T))+log(N)*ncol(B_new);

  WIC_vec_2<-c(WIC_vec_2,WIC_2);
  WIC_vec_log<-c(WIC_vec_log,WIC_log);

  variable.lowest_2<-as.numeric(which(wic1_2==min(wic1_2,na.rm=T))[1]);
  variable.lowest_log<-as.numeric(which(wic1_log==min(wic1_log,na.rm=T))[1]);
  var.low.vec_2<-c(colnames(B_new)[variable.lowest_2+1]);
  var.low.vec_log<-c(colnames(B_new)[variable.lowest_log+1]);
  B_new_2<-as.matrix(B_new[,-(variable.lowest_2+1)]);
  B_new_log<-as.matrix(B_new[,-(variable.lowest_log+1)]);

  cnames_2<-c(cnames_2,list(colnames(B_new_2)));
  cnames_log<-c(cnames_log,list(colnames(B_new_log)));

  for(i in 2:(ncol(B)-1))
    {
    if(i!=(ncol(B)-1))
      {
      wic1_2<-backward_sel_WIC(Y,N,n,B_new_2,id,family,corstr,nb,is.gee);
      wic1_log<-backward_sel_WIC(Y,N,n,B_new_log,id,family,corstr,nb,is.gee);

      wic_mat_2[(i+1),colnames(B_new_2)[-1]]<-wic1_2;
      wic_mat_log[(i+1),colnames(B_new_log)[-1]]<-wic1_log;

      WIC_2<-sum(apply(wic_mat_2[1:(i+1),],1,min,na.rm=T))+2*ncol(B_new_2);
      WIC_log<-sum(apply(wic_mat_log[1:(i+1),],1,min,na.rm=T))+log(N)*ncol(B_new_log);

      WIC_vec_2<-c(WIC_vec_2,WIC_2);
      WIC_vec_log<-c(WIC_vec_log,WIC_log);

      variable.lowest_2<-as.numeric(which(wic1_2==min(wic1_2,na.rm=T))[1]);
      variable.lowest_log<-as.numeric(which(wic1_log==min(wic1_log,na.rm=T))[1]);
      var.low.vec_2<-c(var.low.vec_2,colnames(B_new_2)[variable.lowest_2+1]);
      var.low.vec_log<-c(var.low.vec_log,colnames(B_new_log)[variable.lowest_log+1]);

      B_new_2<-as.matrix(B_new_2[,-(variable.lowest_2+1)]);
      B_new_log<-as.matrix(B_new_log[,-(variable.lowest_log+1)]);
      }

    if(i==(ncol(B)-1))
      {
      if(is.gee==F)
        {
        if(nb==F)
          {
          full.fit_2<-stats::glm(Y~B_new_2-1,family=family);
          full.fit_log<-stats::glm(Y~B_new_log-1,family=family);
          full.wald_2<-(summary(full.fit_2)[12]$coef[-1,3])^2;
          full.wald_log<-(summary(full.fit_log)[12]$coef[-1,3])^2;

          wic1_2<-full.wald_2;
          wic1_log<-full.wald_log;
          }
        if(nb==T)
          {
          full.fit_2<-gamlss::gamlss(Y~B_new_2-1,family="NBI",trace=FALSE);
          full.fit_log<-gamlss::gamlss(Y~B_new_log-1,family="NBI",trace=FALSE);
          sink(tempfile());
          full.wald_2<-((as.matrix(summary(full.fit_2))[,3])[-c(1,nrow(as.matrix(summary(full.fit_2))))])^2;
          full.wald_log<-((as.matrix(summary(full.fit_log))[,3])[-c(1,nrow(as.matrix(summary(full.fit_log))))])^2;
          sink();

          wic1_2<-full.wald_2;
          wic1_log<-full.wald_log;
          }
        }

      if(is.gee==T)
        {
        if(family=="gaussian")
          {
          full.fit_2<-geepack::geeglm(Y~B_new_2-1,id=id,corstr=corstr);
          full.fit_log<-geepack::geeglm(Y~B_new_log-1,id=id,corstr=corstr);

          full.wald_2<-(summary(full.fit_2)[6]$coef[-1,3])^2;
          full.wald_log<-(summary(full.fit_log)[6]$coef[-1,3])^2;
          }
        if(family!="gaussian")
          {
          full.fit_2<-geepack::geeglm(Y~B_new_2-1,id=id,family=family,corstr=corstr);
          full.fit_log<-geepack::geeglm(Y~B_new_log-1,id=id,family=family,corstr=corstr);

          full.wald_2<-(summary(full.fit_2)[6]$coef[-1,3])^2;
          full.wald_log<-(summary(full.fit_log)[6]$coef[-1,3])^2;
          }

        wic1_2<-full.wald_2;
        wic1_log<-full.wald_log;
        }

      wic_mat_2[(i+1),colnames(B_new_2)[-1]]<-wic1_2;
      wic_mat_log[(i+1),colnames(B_new_log)[-1]]<-wic1_log;

      WIC_2<-sum(apply(wic_mat_2[1:(ncol(B)),],1,min,na.rm=T))+2*ncol(B_new_2);
      WIC_log<-sum(apply(wic_mat_log[1:(ncol(B)),],1,min,na.rm=T))+log(N)*ncol(B_new_log);

      WIC_vec_2<-c(WIC_vec_2,WIC_2);
      WIC_vec_log<-c(WIC_vec_log,WIC_log);

      B_new_2<-as.matrix(B_new_2[,-(variable.lowest_2)]);
      B_new_log<-as.matrix(B_new_log[,-(variable.lowest_log)]);

      colnames(B_new_2)<-"Intercept";
      colnames(B_new_log)<-"Intercept";
      }

    cnames_2<-c(cnames_2,list(colnames(B_new_2)));
    cnames_log<-c(cnames_log,list(colnames(B_new_log)));
    }

  graphics::par(mfrow=c(1,2),las=0);
  graphics::plot(rev(WIC_vec_2),type="l",lwd=2,ylab="",xlab="backward path",cex.lab=2);
  graphics::mtext(text="WIC",line=2,side=2,cex=2,las=0);
  graphics::title(main=bquote(paste(~paste(lambda)==.(2))),cex.main=2.2,line=1.4);
  graphics::plot(rev(WIC_vec_log),type="l",lwd=2,ylab="",xlab="backward path",cex.lab=2);
  graphics::mtext(text="WIC",line=2,side=2,cex=2,las=0);
  graphics::title(main=bquote(paste(~paste(lambda)~"= log(N)")),cex.main=2.2,line=1.4);

  if(print.disp==T)
    {
    writeLines("\n Forward pass output: \n");
    forw.info<-cbind(round(score_term,4),pred.name_vec,cut_vec,trunc.type_vec,is.int_vec);
    colnames(forw.info)<-c("Score","Predictor name","Cut term (knot)","No. of new parent terms","Interaction?");
    print(forw.info);
    }

  ## Some final model output, WIC, GCV etc.

  if(is.gee==T)
    {
    B_final<-as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]);
    final_mod_2<-geepack::geeglm(Y~B_final-1,id=id,family=family,corstr=corstr);

    B_final<-as.matrix(B[,colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]);
    final_mod_log<-geepack::geeglm(Y~B_final-1,id=id,family=family,corstr=corstr);
    }
  if(is.gee==F)
    {
    if(nb==F)
      {
      B_final<-as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]);
      final_mod_2<-stats::glm(Y~B_final-1,family=family);
      B_final<-as.matrix(B[,colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]);
      final_mod_log<-stats::glm(Y~B_final-1,family=family);
      }
    if(nb==T)
      {
      B_final<-as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]);
      final_mod.many_2<-mvabund::manyglm(c(t(Y))~B_final-1,family="negative.binomial",maxiter=1000,maxiter2=100);
      final_mod_2<-MASS::glm.nb(c(t(Y))~B_final-1,method="glm.fit2",init.theta=final_mod.many_2$theta);

      B_final<-as.matrix(B[,colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]);
      final_mod.many_log<-mvabund::manyglm(c(t(Y))~B_final-1,family="negative.binomial",maxiter=1000,maxiter2=100);
      final_mod_log<-MASS::glm.nb(c(t(Y))~B_final-1,method="glm.fit2",init.theta=final_mod.many_log$theta);
      }
    }

  NN<-length(Y);

  p_2<-ncol(as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]));
  df1a_2<-p_2+pen*(p_2-1)/2;  # This matches the earth() package, SAS and Friedman (1991) penalty.
  p_log<-ncol(as.matrix(B[,colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]));
  df1a_log<-p_log+pen*(p_log-1)/2;  # This matches the earth() package, SAS and Friedman (1991) penalty.

  RSS1_2<-sum((Y-stats::fitted(final_mod_2))^2);
  RSS1_log<-sum((Y-stats::fitted(final_mod_log))^2);
  RSSq1_2<-1-RSS1_2/sum((Y-mean(Y))^2);
  RSSq1_log<-1-RSS1_log/sum((Y-mean(Y))^2);
  GCV1_2<-RSS1_2/(NN*(1-(df1a_2)/NN)^2);
  GCV1_log<-RSS1_log/(NN*(1-(df1a_log)/NN)^2);

  if(print.disp==T)
    {
    B_final<-as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]);
    final_mod<-final_mod_2;
    wic_mat<-wic_mat_2;
    RSSq1<-RSSq1_2;
    GCV1<-GCV1_2;

    writeLines("\n -- Final model (after pruning/backward selection) for MARGE -- \n");

    if(ncol(B_final)>1){final_mat<-t(t(colnames(as.matrix(B_final))));}
    if(ncol(B_final)==1){final_mat<-as.matrix("Intercept",1,1);}
    colnames(final_mat)<-"Selected variables in the final model:";
    print(final_mat);
    writeLines("\n Final model WIC: \n");
    print(min(wic_mat,na.rm=T));
    writeLines("\n Final model Rsq: \n");
    print(RSSq1);
    writeLines("\n Final model GCV: \n");
    print(GCV1);
    writeLines("\n Final model coefs: \n");
    print(matrix(stats::coef(final_mod),dimnames=list(final_mat)));
    }

  B_final<-list(as.matrix(B[,colnames(B)%in%cnames_2[[which.min(WIC_vec_2)]]]),as.matrix(B[,colnames(B)%in%cnames_log[[which.min(WIC_vec_log)]]]));

  wic_mat<-list(wic_mat_2,wic_mat_log);
  min_wic_own<-list(min(wic_mat_2,na.rm=T),(min(wic_mat_log,na.rm=T)));
  GCV1<-list(GCV1_2,GCV1_log);
  y_pred<-list(stats::predict(final_mod_2),stats::predict(final_mod_log));
  final_mod<-list(final_mod_2,final_mod_log);

  z<-NULL;
  z$bx<-B_final;
  z$wic_mat<-wic_mat;
  z$min_wic_own<-min_wic_own;
  z$GCV<-GCV1;
  z$y_pred<-y_pred;
  z$final_mod<-final_mod;

  class(z)<-"marge";

  return(z);
  }

#' getNumberPart
#'
#' A function to pull out negative and decimals from strings.
#' @name getNumberPart
#' @param x : a numerical value.
#' @export
#' @importFrom gsubfn strapply
getNumberPart<-function(x)
  {
  pat<-"(-?(\\d*\\.*\\d+|\\d+\\.))";
  gsubfn::strapply(x,pattern=pat,FUN=as.numeric,simplify=TRUE,empty=NA);
  }

#' predict.marge
#'
#' A predict-type function for fitted MARS and MARGE objects when using the \code{marge} package.
#' @name predict.marge
#' @param object : the final selected MARS/MARGE model. Only works for \code{mars_ls} and \code{marge} model objects.
#' @param newdata : the new set of predictor variables to make predictions on (test data).
#' @param X_pred : the training predictor variable values from the fitted MARS/MARGE model.
#' @param is.marge : a logical argument, is this a MARGE model? The default is \code{FALSE}.
#' @param pen : the penalty used for the MARGE model (only applicable if \code{marge} was used). The default is \code{pen=2}.
#' @param ... : further arguments passed to or from other methods.
#' @return \code{predict.marge} returns a list of calculated values consisting of:
#' @return \code{eta.p} : the fitted linear predictor using the new (test) data.
#' @return \code{basis_new} : the model matrix for the new (test) data.
#' @author Jakub Stoklosa and David I. Warton.
#' @references Stoklosa, J. and Warton, D.I. (2018). A generalized estimating equation approach to multivariate adaptive regression splines. \emph{Journal of Computational and Graphical Statistics}, \strong{27}, pp. 245--253.
#' @export
#' @importFrom stringr str_extract_all
#' @examples ## Load the "leptrine" presence-absence data.
#'
#' data(leptrine)
#'
#' dat1<-leptrine[[1]]   # Training data.
#' dat1_t<-leptrine[[2]] # Test data.
#'
#' Y<-dat1$Y             # Response variable.
#' N<-length(Y)          # Sample size (number of clusters).
#' n<-1                  # Cluster size.
#' id<-rep(1:N,each=n)   # The ID of each cluster.
#'
#' X_pred<-dat1[,-c(3:10)]    # Design matrix using only two (of nine) predictors.
#' X_predt<-dat1_t[,-c(4:11)]
#'
#' ## Set MARGE tunning parameters.
#'
#' family<-"binomial"    # The selected "exponential" family for the GLM/GEE.
#' is.gee<-FALSE         # Is the model a GEE?
#' nb<-FALSE             # Is this a negative binomial model?
#' tols_score<-0.0001    # A set tolerence (stopping condition) in forward pass for MARGE.
#' M<-21                 # A set threshold for the maximum number of basis functions to be used.
#' print.disp<-FALSE     # Print ALL the output?
#' pen<-2                # Penalty to be used in GCV.
#' minspan<-NULL         # A set minimum span value.
#'
#' ## Fit the MARGE models (about ~ 30 secs.)
#'
#' model_marge<-marge(X_pred,Y,N,n,id,family,corstr,pen,tols_score,M,minspan,print.disp,nb,is.gee)
#'
#' ## Predict on training data.
#'
#' pred_marge_2_y<-predict(model_marge,X_predt,X_pred,TRUE,"2")
#'
#' pred_marge_log_y<-predict(model_marge,X_predt,X_pred,TRUE,"logN")
predict.marge<-function(object,newdata,X_pred,is.marge=F,pen="2",...)
  {
  if(is.marge==T)
    {
    if(pen=="2"){mod<-object$final_mod[[1]];}
    if(pen=="logN"){mod<-object$final_mod[[2]];}
    }

  if(is.marge==F){mod<-object$final_mod;}

  if(ncol(stats::model.matrix(mod))==1){basis_new<-stats::model.matrix(mod);}

  if(ncol(stats::model.matrix(mod))!=1)
    {
    newdat_var<-colnames(newdata)[-1];  # Remove intercept term (1st column).

    fitted_dat<-round(X_pred,digits=4); # Original input data (need the variable names).

    ## Extract the variable names and cuts from the final model basis.

    b1<-substring(colnames(stats::model.matrix(mod)),8)[-1];

    var_name_list<-strsplit(b1,split="\\*");
    qq1<-length(var_name_list);

    pp1<-length(newdat_var);
    pp2<-length(b1);

    basis_new<-rep(1,nrow(newdata));

    for(ww in 1:qq1)
      {
      var1<-var_name_list[[ww]];
      qq2<-length(var1);

      if(qq2==1)  # For additive structures.
        {
        var.vec_new1<-c();

        for(vv in 1:pp1)
          {
          var.vec_new1<-c(var.vec_new1,unique(unlist(stringr::str_extract_all(var1,newdat_var[vv]))));
          }

        var_name<-unlist(stringr::str_extract_all(var1,max(var.vec_new1)));

        lenv1<-length(strsplit(var1,"-|\\s")[[1]]);

        cut1<-c(getNumberPart(var1));

        if(lenv1==2){cut1<-abs(cut1);}
        if(lenv1==3){cut1<-cut1;}

        if(length(cut1)>1 && cut1[1]!=cut1[2])
          {
          if(cut1[1]<0 || cut1[2]<0)
            {
            var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
            cut00<-which(cut1<0);
            cut01<-abs(signif(as.numeric(cut1[cut00]),4));

            cut3<-utils::head(fitted_dat[which(abs(signif(fitted_dat[,var_num0],4))==cut01),var_num0],n=1);

            if(floor(abs(cut1[cut00]))==floor(abs(cut3))){cut2<-abs(cut1[cut00])*sign(cut3);}
            }
          if(cut1[1]>=0 && cut1[2]>=0)
            {
            var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
            cut00<-which(cut1==var_num0);
            cut2=cut3=as.numeric(cut1[-cut00]);
            }
          }
        if(length(cut1)>1 && (abs(cut1[1])==abs(cut1[2])))
          {
          var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
          cut00<-min(cut1)[1];
          cut01<-abs(signif(as.numeric(cut1[cut00]),4));

          cut3<-utils::head(fitted_dat[which(abs(signif(fitted_dat[,var_num0],4))==cut01),var_num0],n=1);
          cut2<-abs(cut1[cut00])*sign(cut3);
          }

        if(length(cut1)==1){cut2=cut3=cut1;}

        temp_name<-paste("(",var_name,"-",cut2,")",sep="");
        var_num<-which(colnames(newdata)==(noquote(var_name)));
        var_chosen<-newdata[,var_num];

        if(temp_name==var1){b1_new<-matrix(tp1(var_chosen,as.numeric(cut3)),ncol=1);}
        if(temp_name!=var1){b1_new<-matrix(tp2(var_chosen,as.numeric(cut3)),ncol=1);}

        basis_new<-cbind(basis_new,round(b1_new,digits=4));
        }

      if(qq2==2)  # For interaction structures.
        {
        basis_new1<-c();

        for(yy in 1:qq2)
          {
          var2<-var1[yy];

          var.vec_new1<-c();

          for(vv in 1:pp1)
            {
            var.vec_new1<-c(var.vec_new1,unique(unlist(stringr::str_extract_all(var2,newdat_var[vv]))));
            }

          var_name<-unlist(stringr::str_extract_all(var2,max(var.vec_new1)));

          lenv1<-length(strsplit(var2,"-|\\s")[[1]]);

          cut1<-c(getNumberPart(var2));

          if(lenv1==2){cut1<-abs(cut1);}
          if(lenv1==3){cut1<-cut1;}

          if(length(cut1)>1 && (cut1[1]!=cut1[2]))
            {
            if(cut1[1]<0 || cut1[2]<0)
              {
              var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
              cut00<-which(cut1<0);
              cut01<-abs(signif(as.numeric(cut1[cut00]),4));

              cut3<-utils::head(fitted_dat[which(abs(signif(fitted_dat[,var_num0],4))==cut01),var_num0],n=1);

              if(floor(abs(cut1[cut00]))==floor(abs(cut3))){cut2<-abs(cut1[cut00])*sign(cut3);}
              }
            if(cut1[1]>=0 && cut1[2]>=0)
              {
              var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
              cut00<-which(cut1==var_num0);
              cut2=cut3=as.numeric(cut1[-cut00]);
              }
            }

          if(length(cut1)>1 && (abs(cut1[1])==abs(cut1[2])))
            {
            var_num0<-which(colnames(fitted_dat)==(noquote(var_name)));
            cut00<-min(cut1)[1];
            cut01<-abs(signif(as.numeric(cut1[cut00]),4));

            cut3<-utils::head(fitted_dat[which(abs(signif(fitted_dat[,var_num0],4))==cut01),var_num0],n=1);
            cut2<-abs(cut1[cut00])*sign(cut3);
            }

          if(length(cut1)==1){cut2=cut3=cut1;}

          temp_name<-paste("(",var_name,"-",cut2,")",sep="");
          var_num<-which(colnames(newdata)==(noquote(var_name)));
          var_chosen<-newdata[,var_num];

          if(temp_name==var2){b1_new<-matrix(tp1(var_chosen,as.numeric(cut3)),ncol=1);}
          if(temp_name!=var2){b1_new<-matrix(tp2(var_chosen,as.numeric(cut3)),ncol=1);}

          basis_new1<-cbind(basis_new1,round(b1_new,digits=4));
          }

        basis_new<-cbind(basis_new,basis_new1[,1]*basis_new1[,2]);
        }
      }
    }

  eta.p<-c(stats::coef(mod)%*%t(basis_new));

  list(eta.p=eta.p,basis_new=basis_new);
  }


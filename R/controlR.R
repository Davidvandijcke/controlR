#-------------------------------------------------------
#  File: controlR.R
#  Purpose: Fully worked-out code for the "Optimal Rearrangement"
#           of the Imbens–Newey control variable, using the 'np' package.
#
#-------------------------------------------------------

#' Estimate the usual Imbens-Newey control variable via nonparametric cCDF
#'
#' This uses the \code{npcdist()} function from the \pkg{np} package
#' to estimate \eqn{F_{X|Z}(x|z)}. Then for each observation \eqn{i},
#' we obtain
#' \deqn{V_i = \hat{F}_{X|Z}( X_i \mid Z_i ).}
#'
#' @param X numeric vector of length n (the endogenous regressor).
#' @param Z numeric matrix/data.frame of instruments (n x d). Multi-dimensional is allowed.
#' @param ... additional options passed to \code{\link[np]{npcdistbw}} or \code{\link[np]{npcdist}}.
#'
#' @return A numeric vector \code{Vhat} of length n (the estimated control variable).
#'
#' @examples
#' \dontrun{
#'   library(np)
#'   set.seed(1)
#'   n <- 200
#'   Z <- matrix(rnorm(n), ncol=1)
#'   X <- Z[,1] + 0.3*rnorm(n)
#'   Vhat <- estimateControlV(X, Z)
#'   hist(Vhat, main="Estimated V = F_{X|Z}(X_i, Z_i)")
#' }
#'
estimateControlV <- function(X, Z, ...) {
  if(!requireNamespace("np", quietly=TRUE)){
    stop("Package 'np' not installed. Please install it to use this function.")
  }
  
  # 1. Bandwidth selection for conditional distribution F_{X|Z}, using npcdistbw
  bwObj <- np::npcdistbw(xdat = Z, ydat = X, ...)
  
  # 2. Fit the conditional distribution object
  cdistObj <- np::npcdist(bws = bwObj, xdat = Z, ydat = X)
  
  # 3. For each i, we want F_{X|Z}( X_i | Z_i )
  #    We do that by predict(..., exdat=Z_i, eydat=X_i).
  #    `exdat` can be a row from Z, `eydat` is the scalar X_i.
  
  n <- length(X)
  Vhat <- numeric(n)
  for(i in seq_len(n)){
    # single row for exdat:
    exz  <- as.data.frame(matrix(Z[i,], nrow=1))
    # single value for X
    eyx  <- X[i]
    # predict cdf
    pred <- np::predict(cdistObj,
                        exdat = exz,
                        eydat = eyx)
    Vhat[i] <- pred
  }
  
  # Vhat should lie in [0,1], but small numerical deviations can happen:
  Vhat <- pmin(pmax(Vhat, 0), 1)
  
  return(Vhat)
}


#' Estimate g_x(u)=f_{V|X}(u|x) and G_x(u)=F_{V|X}(u|x) on a grid of u
#'
#' We fit two nonparametric conditional estimators using \code{np} package:
#' \enumerate{
#'   \item \code{npcdist()} to get \eqn{\widehat{F}_{V|X}(u|x)}
#'   \item \code{npcdens()} to get \eqn{\widehat{f}_{V|X}(u|x)}
#' }
#' Then, for each observation \eqn{i}, and each grid point \eqn{u_k},
#' we evaluate the fitted distribution at \eqn{(x_i, u_k)}.
#'
#' @param X numeric vector (length n).
#' @param V numeric vector (length n), the previously estimated control variable or any \eqn{V}.
#' @param grid_u numeric vector in [0,1] at which we evaluate the cdf/pdf.
#' @param ... extra arguments to pass to the bandwidth selection steps, e.g. \code{npcdistbw}, \code{npcdensbw}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{gx_mat}: (n x length(grid_u)) matrix of \eqn{\hat{f}_{V|X}(u_k | X_i)}.
#'   \item \code{Gx_mat}: (n x length(grid_u)) matrix of \eqn{\hat{F}_{V|X}(u_k | X_i)}.
#' }
#'
#' @examples
#' \dontrun{
#'   library(np)
#'   n <- 200
#'   # Suppose we have X and V from somewhere
#'   X <- rnorm(n)
#'   V <- pnorm(X + 0.2*rnorm(n)) # toy
#'   grid_u <- seq(0,1,length.out=21)
#'   est <- estimate_gx_Gx(X, V, grid_u)
#'   dim(est$gx_mat)
#'   dim(est$Gx_mat)
#' }
estimate_gx_Gx <- function(X, V, grid_u=seq(0,1,length.out=50), ...){
  if(!requireNamespace("np", quietly=TRUE)){
    stop("Package 'np' not installed.")
  }
  n <- length(X)
  
  # 1. Fit cdf for V|X
  #    We do bandwidth selection + model fit.
  #    Then we can evaluate F_{V|X}(u|x).
  
  # Note: for multi-dim X, pass a matrix X. 
  #       here X is numeric vector, so that is 1D or unidim with shape (n,1).
  cdist_bw <- np::npcdistbw(xdat=X, ydat=V, ...)
  cdistObj <- np::npcdist(bws=cdist_bw, xdat=X, ydat=V)
  
  # 2. Fit pdf for V|X
  #    We'll use npcdens, which does conditional density f_{V|X}(v|x).
  cdens_bw <- np::npcdensbw(xdat=X, ydat=V, ...)
  cdensObj <- np::npcdens(bws=cdens_bw, xdat=X, ydat=V)
  
  # 3. For each i and each grid_u[k], evaluate cdf & pdf
  #    We'll store in Gx_mat[i,k], gx_mat[i,k].
  
  m <- length(grid_u)
  Gx_mat <- matrix(NA_real_, nrow=n, ncol=m)
  gx_mat <- matrix(NA_real_, nrow=n, ncol=m)
  
  for(i in seq_len(n)){
    # make data.frame for exdat
    # if X is vector, we do a single row with X_i
    Xi_df <- as.data.frame(matrix(X[i], nrow=1))
    
    for(k in seq_len(m)){
      u_k <- grid_u[k]
      # Evaluate cdf:  F_{V|X}(u_k | X_i)
      #  predict( cdistObj, exdat= Xi_df, ey dat = u_k )
      cdf_val <- np::predict(cdistObj,
                             exdat = Xi_df,
                             eydat = u_k)
      
      # Evaluate pdf:  f_{V|X}(u_k | X_i)
      # for npcdens, we pass newdata with columns (X, V).
      newdata_df <- data.frame(X = X[i], V = u_k)
      pdf_val <- np::predict(cdensObj, newdata = newdata_df)
      
      Gx_mat[i,k] <- cdf_val
      gx_mat[i,k] <- pdf_val
    }
  }
  
  return(list(gx_mat=gx_mat, Gx_mat=Gx_mat))
}


#' Compute Q^*_R(u) and do the monotonic rearrangement
#'
#' We form
#' \deqn{f_V(u) = \frac{1}{n}\sum_i g_{x_i}(u),\quad
#'       Q^*_R(u) = \frac{1}{n}\sum_i \left[\frac{g_{x_i}(u)}{f_V(u)}\cdot G_{x_i}(u)\right].}
#'
#' Then we perform an isotonic regression on \eqn{u\mapsto Q^*_R(u)} to get the
#' non-decreasing rearrangement \eqn{\widetilde{Q}_R^*(u)}.
#'
#' @param gx_mat matrix of dimension n x \code{length(grid_u)}, for \eqn{g_{x_i}(u_k)}
#' @param Gx_mat matrix of dimension n x \code{length(grid_u)}, for \eqn{G_{x_i}(u_k)}
#' @param grid_u numeric vector in [0,1] giving the evaluation points
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{Qraw}: numeric vector of length(grid_u) = the unconstrained \eqn{Q^*_R(u_k)}
#'   \item \code{Qmono}: numeric vector of length(grid_u) = the rearranged monotone version
#'   \item \code{grid_u}: the same grid
#' }
#'
#' @examples
#' # Suppose we already have gx_mat, Gx_mat from estimate_gx_Gx().
#' # We'll just show the shape:
#' \dontrun{
#'    qres <- computeQRstar(gx_mat, Gx_mat, grid_u)
#'    plot(grid_u, qres$Qraw, type='l', main="Q^*_R(u)")
#'    lines(grid_u, qres$Qmono, col=2, lty=2, lwd=2)
#' }
computeQRstar <- function(gx_mat, Gx_mat, grid_u){
  n <- nrow(gx_mat)
  m <- length(grid_u)
  
  # unconditional pdf of V: f_V(u) = (1/n) * sum_i gx_mat[i,k]
  fV_vec <- colMeans(gx_mat)
  
  Qraw <- numeric(m)
  for(k in seq_len(m)){
    denom <- fV_vec[k]
    if(denom < 1e-16) {
      Qraw[k] <- 0
    } else {
      ratio_i <- gx_mat[,k]/denom
      Qraw[k] <- mean( ratio_i * Gx_mat[,k] )
    }
  }
  
  # monotonic rearrangement using isotonic regression
  isoRes <- isoreg(x=grid_u, y=Qraw)
  Qmono <- isoRes$yf
  
  return(list(Qraw=Qraw, Qmono=Qmono, grid_u=grid_u))
}


#' Main function: Construct the "Optimal Rearranged Control Variable" via np
#'
#' This function:
#' \enumerate{
#'   \item Estimates V_i = F_{X|Z}(X_i,Z_i) using \code{estimateControlV}.
#'   \item Estimates g_x(u) and G_x(u) at grid_u using \code{estimate_gx_Gx}.
#'   \item Computes Q^*_R(u) and does monotonic rearrangement (using \code{computeQRstar}).
#'   \item Maps each V_i -> R_i by interpolating \eqn{\widetilde{Q}_R^*} at V_i.
#' }
#'
#' @param X numeric vector (length n), the endogenous regressor.
#' @param Z numeric matrix/data.frame (n x d) of instruments.
#' @param grid_u numeric vector of points in [0,1] for rearrangement (default length=50).
#' @param ... additional arguments passed to the bandwidth selection in each sub-step,
#'        e.g. \code{bwmethod="cv.ls"} or others from \pkg{np}.
#'
#' @return A list:
#' \itemize{
#'   \item \code{V}: the usual Imbens--Newey control variable estimates
#'   \item \code{R}: the optimally rearranged control variable (length n)
#'   \item \code{Qraw}, \code{Qmono}: the raw and rearranged \eqn{Q_R(u)} on the grid
#'   \item \code{grid_u}: the grid used
#' }
#'
#' @examples
#' \dontrun{
#'   library(np)
#'   set.seed(123)
#'   n <- 100
#'   Z <- rnorm(n)
#'   X <- Z + 0.5*rnorm(n)   # Triangular data generating process
#'
#'   out <- controlR_np(X, Z, grid_u=seq(0,1,length.out=40))
#'   # Now out$R is the rearranged control. Compare hist:
#'   par(mfrow=c(1,2))
#'   hist(out$V, main="Usual I-N Control V", xlim=c(0,1))
#'   hist(out$R, main="Optimal Rearranged Control R", xlim=c(0,1))
#'
#'   # Next step: run a nonparametric regression of Y on (X,R).
#' }
controlR <- function(X, Z, grid_u=seq(0,1,length.out=50), ...){
  # 1. Estimate the usual Imbens-Newey control
  Vhat <- estimateControlV(X, Z, ...)
  
  # 2. Estimate g_x_i(u), G_x_i(u) on a grid of u
  est <- estimate_gx_Gx(X, Vhat, grid_u=grid_u, ...)
  gx_mat <- est$gx_mat
  Gx_mat <- est$Gx_mat
  
  # 3. Compute Q^*_R(u) and do monotonic rearrangement
  qres <- computeQRstar(gx_mat, Gx_mat, grid_u)
  
  # 4. Interpolate R_i = Qmono( V_i ).  We'll do simple linear interpolation.
  Qmono <- qres$Qmono
  Rhat <- numeric(length(X))
  
  for(i in seq_along(X)){
    vi <- Vhat[i]
    if(vi <= grid_u[1]){
      Rhat[i] <- Qmono[1]
    } else if(vi >= grid_u[length(grid_u)]){
      Rhat[i] <- Qmono[length(grid_u)]
    } else {
      idx <- findInterval(vi, grid_u)
      # linear interpolation between (grid_u[idx], Qmono[idx]) and (grid_u[idx+1], Qmono[idx+1])
      x0 <- grid_u[idx]
      x1 <- grid_u[idx+1]
      y0 <- Qmono[idx]
      y1 <- Qmono[idx+1]
      w  <- (vi - x0)/(x1 - x0)
      Rhat[i] <- (1-w)*y0 + w*y1
    }
  }
  
  list(
    V      = Vhat,            # usual control
    R      = Rhat,            # rearranged control
    Qraw   = qres$Qraw,       # raw Q_R^*(u)
    Qmono  = Qmono,           # rearranged monotone
    grid_u = qres$grid_u
  )
}

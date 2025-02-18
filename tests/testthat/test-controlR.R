###########################################################
# Test code with simulated data for "optimal rearrangement"
###########################################################

# 1. Load or source your "optimalControl.R" file:
#    Make sure 'controlR_np()', 'estimateControlV_np()', etc. are accessible.
devtools::load_all()  # Adjust path if needed.

# 2. (Optional) Load the np package explicitly
#    install.packages("np")  # if not installed
library(np)

#-----------------------------------------------------------
# 3. Simulate data from a simple triangular model
#    We'll do:
#       Z ~ N(0,1)
#       eta ~ N(0,1)
#       X = Z + 0.5*eta  (strictly increasing in eta, so valid triangular)
#       eps ~ N(0,1)
#       Y  = 2 * sin(X) + eps   (an arbitrary nonparametric function)
# 
#    This structure is purely illustrative.
#-----------------------------------------------------------
set.seed(123)
n <- 2000
Z <- rnorm(n)
eta <- rnorm(n)
X <- Z + 0.5 * eta   # "Reduced form" with endogeneity
eps <- rnorm(n)
Y <- 2 * sin(X) + eps   # Some structural equation for the outcome

#-----------------------------------------------------------
# 4. Fit the "optimal rearranged" control
#    - This uses controlR_np(), which in turn uses:
#         estimateControlV_np()  (for V = F_{X|Z}(X,Z))
#         estimate_gx_Gx_np()    (for g_x, G_x)
#         computeQRstar_np()     (for Q_R*(u) + rearrangement)
#-----------------------------------------------------------
out <- controlR(X, Z, grid_u=seq(0,1, length.out=40),
                   # pass optional arguments for bandwidth selection:
                   bwmethod = "cv.ls",  # e.g., cross-validation
                   ckertype = "epanechnikov" 
)

# test
# out$V = usual Imbens–Newey control
# out$R = optimally rearranged control
cat("Summary of V:\n")
print(summary(out$V))
cat("\nSummary of R:\n")
print(summary(out$R))

# Check how V and R differ in distribution:
op <- par(mfrow=c(1,2))
hist(out$V, breaks=15, main="Histogram: V", xlab="V", col="lightblue")
hist(out$R, breaks=15, main="Histogram: R (Rearranged)", xlab="R", col="lightgreen")
par(op)

# Possibly we want to see Qraw(u) vs Qmono(u):
plot(out$grid_u, out$Qraw, type='l',
     main="Q_R*(u) raw vs. monotonic rearrangement",
     xlab="u in [0,1]", ylab="Q_R*(u)")
lines(out$grid_u, out$Qmono, col=2, lty=2, lwd=2)
legend("topleft", legend=c("Qraw","Qmono"), col=c(1,2), lty=c(1,2))

#-----------------------------------------------------------
# 5. Illustrate second-stage usage:
#    We can estimate E[Y | X, R] by a standard nonparametric regression
#    using (X, out$R) as regressors. We'll show a quick example with 'npreg'.
#-----------------------------------------------------------
fit2 <- npreg(bws = NULL, # let np auto-select or specify
              tydat = Y,
              xdat = data.frame(X=X, R=out$R), 
              # optionally pass e.g. ckertype="epanechnikov", etc.
)
summary(fit2)

# We can do some quick predictions on a grid:
xseq <- seq(min(X), max(X), length.out=50)
# We'll just fix R=median(out$R) for illustration
rfix <- median(out$R)
pred_grid <- npreg:::predict.npreg(fit2,
                                   newdata = data.frame(X=xseq, R=rep(rfix,50)))
plot(xseq, pred_grid, type='l',
     main="Estimated E[Y | X=x, R=median(R)]",
     xlab="X", ylab="predicted Y")

cat("\nDone.\n")

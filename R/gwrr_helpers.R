# Helper functions for GWRR

# Do bi-section section search for phi
gwr.bw.cv <- function(lb, ub, eps, X, y, S, kernel){
   a <- lb
   b <- ub
   c <- (a+b)/2
   diff <- b - a
   N <- dim(X)[1]

   while (diff > eps){
      a.c <- (a+c)/2
      c.b <- (c+b)/2
      RMSE.a.c <- gwr.cv.err(a.c, X, y, S, N, kernel)
      RMSE.c.b <- gwr.cv.err(c.b, X, y, S, N, kernel)

      if (RMSE.a.c < RMSE.c.b){
         b <- c.b
         RMSE.b <- RMSE.c.b
         print(paste("Bandwidth: ", format(b,digits=4), " RMSPE :", format(RMSE.b,digits=4)))
      }

      if (RMSE.a.c > RMSE.c.b){
         a <- a.c
         RMSE.a <- RMSE.a.c
         print(paste("Bandwidth: ", format(a,digits=4), " RMSPE :", format(RMSE.a,digits=4)))
      }

      c <- (a+b)/2
      diff <- abs(b - a)
   }

   RMSE.lb <- gwr.cv.err(lb, X, y, S, N, kernel)
   RMSE.ub <- gwr.cv.err(ub, X, y, S, N, kernel)
   RMSE.c <- gwr.cv.err(c, X, y, S, N, kernel)

   # Check bounds
   if (RMSE.lb < RMSE.c){
      c <- lb
      RMSE.c <- RMSE.lb
   }
   if (RMSE.ub < RMSE.c){
      c <- ub
      RMSE.c <- RMSE.ub
   }
   print(paste("Bandwidth: ", format(c,digits=4), " RMSPE :", format(RMSE.c,digits=4)))
   params <- list(c, RMSE.c, RMSE.c * N)
   names(params) <- c("phi", "RMSPE", "cv.score")
   params
}

# Calculate CV RMSPE for all observations
gwr.cv.err <- function(phi, X, y, S, N, kernel){
   yhat <- array(0,N)
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)

   for (i in 1:N){
      W.i <- diag(W[i,])
      W.i[i,i] <- 0
      beta <- chol2inv(chol(t(X) %*% W.i %*% X)) %*% t(X) %*% W.i %*% y
      yhat[i] <- X[i,] %*% beta
   }
   sqrt(sum((y - as.vector(yhat))^2) / N)
}

# Function to calculate spatial correlation
w.exp <- function(phi, S) exp(-(S/phi))

# Function to calculate spatial correlation
w.gauss <- function(phi, S) exp(-(S/phi)^2)

# Estimate GWR beta coefficients given bandwidth
gwr.beta <- function(phi, X, y, S, N, kernel){
   if (kernel == "exp") W <- w.exp(phi, S)
   if (kernel == "gauss") W <- w.gauss(phi,S)

   p <- dim(X)[2]
   betas <- array(0, dim=c(p,N))

   # Estimate regression coefficients
   for (i in 1:N){
      W.i <- diag(W[i,])
      betas[,i] <- t(chol2inv(chol(t(X) %*% W.i %*% X)) %*% t(X) %*% W.i %*% y)
   }
   betas
}

gwr.yhat <- function(betas, X){
   N <- dim(X)[1]
   yhat <- array(0, N)
   for (i in 1:N) yhat[i] <- X[i,] %*% betas[,i]
   yhat
}

gwr.rmse <- function(y, yhat){
   N <- length(y)
   sqrt(sum((y - yhat)^2) / N)
}

gwr.rsquare <- function(y, yhat){
   rs1 <- (cor(y, yhat))^2
   SS.err <- sum((y - yhat)^2)
   ybar <- mean(y)
   SS.tot <- sum((y - ybar)^2)
   rs2 <- 1 - (SS.err/SS.tot)
   rs2
}

# Calculate the variance decomposition proportions and scaled condition indexes at one location
local.vdp <- function(W.X, N){
   # Scale input matrix W.X by norm first
   one <- rep(1,N)
   normx <- sqrt(drop(one %*% (W.X^2)))
   W.X <- scale(W.X, FALSE, normx)

   # svd returns: P is vector of singular values,
   # U is matrix of left singular vectors, V is matrix of right singular vectors
   svd.list <- svd(W.X)
   P <- svd.list$d
   #U <- svd.list$u
   V <- svd.list$v
   p <- length(P)
   phi.ij <- array(0, dim = c(p,p))
   phi.i <- array(0,p)
   pi.ij <- array(0, dim = c(p,p))
   k.X <- array(0, dim=c(p,1))

   #Use the X'X = VS^2V' result
   for (j in 1:p){
      phi.ij[,j] <- V[,j]^2 / P[j]^2
   }

   phi.i <- apply(phi.ij, 1, sum)

   for (i in 1:p){
      pi.ij[,i] <- phi.ij[i,] / phi.i[i]
   }

   for (i in 1:p){
      k.X[i] <- P[1] / P[i]
   }

   params <- list(k.X, pi.ij)
   names(params) <- c("k.X","pi.ij")
   params
}

#'estimate covariance function
#'@export
#'@param y  n*p time series
estCov <- function(y, k0 = 3)
{
  Covfun <- function(ym, k = 1)
  {
    n <- nrow(ym)
    p <- ncol(ym)
    Covk <- matrix(0, p, p)
    ybar <- colMeans(ym)
    for(i in 1:(n-k))
    {
      Covk <- Covk + (ym[i+k,] - ybar) %*% t(ym[i,] - ybar)
    }

    Covk <- Covk/(n-k)
    return(Covk)
  }



  Mhat <- matrix(0, ncol(y), ncol(y))

  for(j in 1:k0)
  {
    Sj <- Covfun(y, k = j)
    Mhat <- Mhat + Sj %*% t(Sj)
  }

  return(Mhat)

}


estCorr <- function(y, k0 = 3)
{
  Covfun <- function(ym, k = 1)
  {
    n <- nrow(ym)
    p <- ncol(ym)
    Covk <- matrix(0, p, p)
    ybar <- colMeans(ym)
    for(i in 1:(n-k))
    {
      Covk <- Covk + (ym[i+k,] - ybar) %*% t(ym[i,] - ybar)
    }

    Covk <- Covk/(n-k)
    return(Covk)
  }



  Mhat <- matrix(0, ncol(y), ncol(y))

  for(j in 1:k0)
  {
    Sj <- Covfun(y, k = j)
    Mhat <- Mhat + Sj %*% t(Sj)
  }

  vary <- diag(Mhat)

  Mhat <- Mhat/(sqrt(vary)) %*% t(sqrt(vary))

  return(Mhat)

}



#'a function to find a
#'@export
#'@param lambda eigenvalues

findr <- function(lambda)
{
  p <- length(lambda)
  R <- p/2
  rhat <- which.min(lambda[2:(R+1)]/lambda[1:R])
  return(rhat)
}

#mcp
mcp <- function(x, gam = 3, lambda)
{
  if(abs(x) <= gam*lambda)
  {
    value <- lambda*abs(x) - x^2/(2*gam)
  }
  else{
    value <- 0.5*gam*lambda^2
  }
}
mcp <- Vectorize(mcp,"x")

#mcp vector
mcpv <- function(x, gam = 3, lambda)
{
  abs_x <-  sqrt(sum(x^2))

  if(abs_x <= gam*lambda)
  {
    value <- lambda*abs_x - abs_x^2/(2*gam)
  }else{
    value <- gam*lambda^2/2
  }
  return(value)
}


#sfunction
sfun <- function(x, th)
{
  xn <- sqrt(sum(x^2))

  if(xn == 0)
  {
    return(x)
  }
  thval <- 1 - th/xn
  #thval*((thval) >0)*x
  max(thval, 0)*x
}

#mcp solution
mcp_sol <- function(x,lam,nu,gam)
{
  temp <- gam*lam
  xn <- sqrt(sum(x^2))
  if(xn <= temp)
  {
    z <- sfun(x,lam/nu) / (1 - 1/(gam*nu))
  }else{
    z <- x
  }
  return(z)
}



#'fitted value
#'@export
fittedfun <- function(Thetaest,ym)
{
  n <- nrow(ym)
  p <- ncol(ym)
  predm <- matrix(0, n, p)

  for(i in 1:n)
  {
    fi <- solve(t(Thetaest) %*% Thetaest) %*% t(Thetaest) %*% ym[i,]
    predm[i,] <-  Thetaest %*% fi
  }
  return(predm)
}


#'get estimated factor process
#'@export
#'@param obj estimated object
#'@param ym observed time series
estfactor <- function(obj,ym)
{
  n <- nrow(ym)
  Thetaest <- obj$fit$Qest
  fm <- matrix(0, n, ncol(Thetaest))
  for(i in 1:n)
  {
    fm[i,] <- solve(t(Thetaest) %*% Thetaest) %*% t(Thetaest) %*% ym[i,]
  }
  return(fm)
}



#'  Sparse factor models
#' @import MASS
#' @import ncvreg
#' @import Matrix
#' @param Shat the traditional estimator: a p*r matrix
#' @param r number of factors
#' @param lambda a given lambda value
#' @return
#' @export

sparsefa_seq_Q <- function(Qt, r, penalty ="MCP", lambda, rho = 1,
                  niters = 1000, tol1 = 1e-6, seed = 123)
{
  set.seed(seed)
  p <- nrow(Qt)
  A <- Qt %*% t(Qt)
  B <- diag(p)

  Ip <- diag(1,p)

  Sm <- matrix(0, p, r)
  Qest <- matrix(0,p, r)
  objvalue <- rep(0,r)

  ### initial values
  v <- rep(0, p)
  d1_old <- 1

  if(r>=2)
  {
    varimaxQ <- unclass(varimax(Qt)$loadings)
    varimaxQ <- varimaxQ[,sort(colSums(abs(varimaxQ)), index.return = TRUE)$ix]
    q0 <- varimaxQ[,1]
    q0 <- ginv(B) %*% q0
  }else{
    q0 <- ginv(B) %*% Qt[,1]
  }


  ##### first component
  resj <- seqfunQ(A, B, p,q0,d1_old, v, niters, lambda, tol1, rho, penalty)

  Sm[,1] <- resj$zm
  Qest[,1] <- resj$qm * sign(resj$qm[which(round(resj$qm,5)!=0)[1]])

  obj <- 0.5*norm(A - Sm %*% t(Sm),"F")^2 +
    sum(unlist(lapply(Qest[,1],FUN = function(x){mcp(x,3,lambda)})))
  objvalue[1] <- obj


  ### rth components
  if(r >=2)
  {
    v <- rep(0, p)

    for(j in 2:r)
    {
      B <- as.matrix(forceSymmetric(Ip - Sm %*% t(Sm))) ## new B

      q0 <- ginv(B) %*% varimaxQ[,j]

      resj <- seqfunQ(A, B, p,q0, d1_old, v, niters, lambda, tol1, rho, penalty)

      Sm[,j] <- resj$zm
      Qest[,j] <- resj$qm


      obj <- 0.5*norm(A - resj$zm %*% t(resj$zm),"F")^2 +
        sum(unlist(lapply(Qest[,j],FUN = function(x){mcp(x,3,lambda)})))
      objvalue[j] <-  obj
    }
  }

  obj0 <- 0.5*norm(A - Sm %*% t(Sm),"F")^2 +
    sum(unlist(lapply(Qest,FUN = function(x){mcp(x,3,lambda)})))

  Qest <- round(Qest, 5)

  return(list(Sest = Sm,  Qest = Qest,
              objvalue = objvalue, obj0 = obj0,
              niters = resj$niters, lambda = lambda))

}




seqfunQ <- function(A, B, p, q_old, d1_old, v, niters, lam1, tol1, rho, penalty)
{
  for(mm in 1:niters)
  {
    ## update z

    C2 <- A %*% B %*% q_old * d1_old + B%*% q_old * rho -v
    z_updated <- C2/norm(C2,"2")

    # update q

    # yq <- (v + rho * z_updated + d1_old * A %*% z_updated)/sqrt(((rho + d1_old^2)))
    # Xq <- as.matrix(sqrt(((rho + d1_old^2))) * B)

    yq <- sqrt(p)*(v + rho * z_updated + d1_old * A %*% z_updated)/sqrt(((rho + d1_old^2)))
    Xq <- sqrt(p)*as.matrix(sqrt(((rho + d1_old^2))) * B)

    q_updated <- ncvfit(X = Xq , y = yq,init = q_old, penalty = penalty,
                        lambda = lam1, warn = FALSE)$beta

    v_updated <- v + rho * (z_updated - B %*% q_updated)

    diff_norm = norm(z_updated - B %*% q_updated,"2")

    z_old <- z_updated
    q_old <- q_updated
    v <- v_updated

    if(diff_norm <= tol1)
    {
      break
    }

    #rho = rho
  }

  zm <- B%*%q_updated
  zm <- zm/norm(zm,"2")

  return(list(zm = zm, qm = q_updated, niters = mm))
}

#'  Sparse factor models using bic
#' @export
#' @param ym a n*p time series matrix
#' @param r number of factors
#' @param lambda a vector of pre-defined lambda value
#'
sparsefa_seq_bic <- function(ym, r, lambda)
{
  p <- ncol(ym)
  n <- nrow(ym)

  Mhat <- estCov(ym, 1)
  Qt <- as.matrix(eigen(Mhat)$vectors[,1:r])

  nlam1 <- length(lambda)
  bic1 <- rep(0, nlam1)
  loss1 <- rep(0, nlam1)
  param1 <- rep(0, nlam1)

  for(i in 1:nlam1)
  {
    resi <- sparsefa_seq_Q(Qt = Qt,r = r,lambda = lambda[i])
    if(mean(round(abs(resi$Qest),5)>0) >0)
    {
      predi <- fittedfun(resi$Qest,ym)
      loss1[i] <- mean((ym - predi)^2)
      param1[i] <- sum(round(resi$Qest,5)!=0)
      bic1[i] <- log(loss1[i]) + param1[i]*log(n*p)/n/p
    }
    else{break}
  }


  res1 <- sparsefa_seq_Q(Qt = Qt,r = r,lambda = lambda[which.min(bic1)])

  out <- list(fit = res1, bic = bic1, lambda = lambda, lammin = lambda[which.min(bic1)])
  return(out)

}





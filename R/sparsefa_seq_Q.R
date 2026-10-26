#'  Sparse factor models
#' @import MASS
#' @import ncvreg
#' @import Matrix
#' @param Qt the traditional estimator: a p*r matrix
#' @param r number of factors
#' @param lambda a given lambda value
#' @return A list with components: \code{Sest}, the estimated p*r loading
#'   matrix; \code{Qest}, the estimated sparse loading matrix with
#'   all-zero columns dropped; \code{objvalue}, a length-r vector of objective
#'   values per component; \code{obj0}, the overall objective value;
#'   \code{niters}, the number of iterations used by the final component; and
#'   \code{lambda}, the tuning parameter used.
#' @export

sparsefa_seq_Q <- function(Qt, r, penalty ="MCP", lambda, varimaxQ,
                           rho = 1, niters = 1000, tol1 = 1e-6, seed = 123)
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
    q0 <- varimaxQ[,1]
  }else{
    q0 <- Qt[,1]
  }


  ##### first component
  resj <- seqfunQ(A, B, p, q0, d1_old, v, niters, lambda, tol1, rho, penalty)

  Sm[,1] <- resj$zm
  Qest[,1] <- resj$qm

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

  Qest <- Qest[,colSums(abs(Qest)==0)<p,drop = FALSE]

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
  }

  if(norm(q_updated, "2")!=0)
  {
    normz <- norm(B%*%q_updated,"2")
    q_updated <- q_updated/normz
  }

  zm <- B %*% q_updated

  return(list(zm = zm, qm = q_updated, niters = mm))
}


#'  Sparse factor models using bic
#' @export
#' @param ym a n*p time series matrix
#' @param r number of factors
#' @param lambda a vector of pre-defined lambda value
#' @param penalty penalty option in ncvreg, can be MCP, SCAD, lasso
#'
sparsefa_seq_bic <- function(ym, r, lambda, h0=1, penalty = "MCP")
{
  p <- ncol(ym)
  n <- nrow(ym)

  Mhat <- estCov(ym, h0)
  Qt <- as.matrix(eigen(Mhat)$vectors[,1:r])

  nlam1 <- length(lambda)
  bic1 <- rep(NA, nlam1)
  loss1 <- rep(0, nlam1)
  param1 <- rep(0, nlam1)

  varimaxQ <- unclass(varimax(Qt)$loadings)
  varimaxQ <- varimaxQ[,sort(colSums(abs(varimaxQ)), index.return = TRUE)$ix]

  for(i in 1:nlam1)
  {
    resi <- sparsefa_seq_Q(Qt = Qt,r = r,lambda = lambda[i], penalty = penalty,
                           varimaxQ = varimaxQ)
    if(ncol(resi$Qest) == 0){break}
    predi <- fittedfun(resi$Qest,ym)
    loss1[i] <- mean((ym - predi)^2)
    param1[i] <- sum(round(resi$Qest,10)!=0)
    bic1[i] <- log(loss1[i]) + param1[i]*log(n*p)/n/p
  }


  res1 <- sparsefa_seq_Q(Qt = Qt,r = r,lambda = lambda[which.min(bic1)],
                         penalty = penalty, varimaxQ = varimaxQ)


  out <- list(fit = res1, bic = bic1, lambda = lambda,
              lammin = lambda[which.min(bic1)])
  return(out)
}





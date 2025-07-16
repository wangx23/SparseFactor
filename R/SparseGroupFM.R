#'  Sparse Group Factor models
#' @import MASS
#' @import ncvreg
#' @import Matrix
#' @return
#' @export
#' @param Shat the traditional estimator: a p*r matrix
#' @param r number of factors
#' @param group  a p*r matrix, to indicate the group structure for each loading


SparseGroupFM <- function(Shat, r, group, penalty ="MCP",
                          lambda1, lambda2, rho1 = 1, rho2 = 1, gam = 3,
                          niters = 1000, tol = 1e-6)
{
  p <- nrow(Shat)
  G <- Shat %*% t(Shat)
  B <- diag(p)

  Ip <- diag(1,p)

  Sest <- matrix(0, p, r) ## estimated S
  Qest <- matrix(0,p, r) ### estimated Q
  objvalue <- rep(0,r)

  ### initial values
  ### q0 for first component ###
  if(r>=2)
  {
    varimaxS <- unclass(varimax(Shat)$loadings)
    varimaxS <- varimaxS[,sort(colSums(abs(varimaxS)), index.return = TRUE)$ix]
    q0 <- varimaxS[,1]
  }else{
    q0 <- Shat[,1]
  }

  q0 <- q0*sign(q0[1]) ### first component positive

  ##### first component
  group_vec1 <- group[,1]
  fit1 <- seqfunQ_sparse(G, B, p, q0, group_vec1,
                  lam1 =lambda1, lam2 = lambda2,
                  rho1, rho2, gam, penalty, niters,tol)


  Sest[,1] <- fit1$sm
  Qest[,1] <- fit1$qm

  obj <- 0.5*norm(G - Sest %*% t(Sest),"F")^2 +
    sum(unlist(lapply(Qest[,1],FUN = function(x){mcp(x,gam,lambda1)}))) +
    sum(unlist(lapply(1:max(group_vec1), function(x){mcpv(Qest[group_vec1==x,1], gam,
                                                          sqrt(sum(group_vec1==x))*lambda2)})))

  objvalue[1] <- obj


  ### rth components
  if(r >=2)
  {
    for(j in 2:r)
    {
      Bj <- as.matrix(forceSymmetric(Ip - Sest %*% t(Sest))) ## new B
      q0j <- ginv(Bj) %*% varimaxS[,j]
      group_vecj <- group[,j]

      fitj <- seqfunQ_sparse(G, Bj, p, q0j, group_vecj,
                      lam1 =lambda1, lam2 = lambda2,
                      rho1, rho2, gam, penalty, niters, tol)


      Sest[,j] <- fitj$sm
      Qest[,j] <- fitj$qm


      obj <- 0.5*norm(G - fitj$sm %*% t(fitj$sm),"F")^2 +
        sum(unlist(lapply(Qest[,j],FUN = function(x){mcp(x,gam,lambda1)}))) +
        sum(unlist(lapply(1:max(group_vecj), function(x){mcpv(Qest[group_vecj==x,j], gam,
                                                              sqrt(sum(group_vecj==x))*lambda2)})))

      objvalue[j] <-  obj
    }
  }

  obj_p1 <- 0.5*norm(G - Sest %*% t(Sest),"F")^2
  obj_p2 <- sum(unlist(lapply(Qest,FUN = function(x){mcp(x,gam, lambda1)})))
  obj_p3 <- 0

  for(j in 1:r)
  {
    obj_p3 <- obj_p3 +
      sum(unlist(lapply(1:max(group[,j]), function(x){mcpv(Qest[group[,j]==x,j],
                                                           gam, sqrt(sum(group[,j]==x))*lambda2)})))
  }
  obj0 <-  obj_p1 + obj_p2 + obj_p3

  Qest <- round(Qest, 5)

  return(list(Sest = Sest,  Qest = Qest,
              objvalue = objvalue, obj0 = obj0,
              niters = fitj$niters,
              lambda1 = lambda1, lambda2= lambda2))

}




seqfunQ_sparse <- function(G, B, p, q0, group_vec, lam1, lam2,
                    rho1, rho2, gam, penalty, niters, tol)
{
  J <- max(group_vec)
  nJ <- rep(0,J)
  for(j in 1:J)
  {
    nJ[j] <- sum(group_vec==j)
  }

  ### initial values
  delta0 <- q0
  v1 <- v2 <- rep(0,p)

  q_cur <- q0
  delta_cur <- delta0

  Ip <- diag(1, p)
  Xq <- chol(rho2* B + rho1*Ip)
  yq1 <- solve(t(Xq))

  delta_update <- delta_cur

  for(mm in 1:niters)
  {
    ## update s
    C2 <- G %*% B %*% q_cur + rho1* B %*% q_cur -v1
    s_update <- C2/norm(C2,"2")

    # update q

    yq <- yq1 %*% (B %*% G %*% s_update + rho1 * B%*% s_update + B%*% v1 + rho2*delta_cur + v2)
    q_update <- ncvfit(X = sqrt(p)*Xq , y = sqrt(p)*yq,init = q_cur, penalty = penalty,
                        lambda = lam1, warn = FALSE)$beta

    # update delta
    uval <- q_update - v2/rho2

    for(j in 1:J)
    {
      indexj <- group_vec ==j
      delta_update[indexj] <- mcp_sol(uval[indexj],sqrt(nJ[j])*lam2, rho2, gam)
    }

    # update v
    v1 <- v1 + rho1 * (s_update - B %*% q_update)
    v2 <- v2 + rho2 * (delta_update - q_update)

    diff_norm1 <- norm(s_update - B %*% q_update,"2")
    diff_norm2 <- norm(delta_update - q_update, "2")

    q_cur <- q_update
    delta_cur <- delta_update

    if(diff_norm1 + diff_norm2 <= sqrt(p)*tol)
    {
      break
    }

  }

  sm <- B %*% q_update
  sm <- sm/norm(sm,"2")

  return(list(sm = sm, qm = q_update, niters = mm))
}




#'  Sparse Group Factor models based on BIC, same lambda for two tuning
#' @import MASS
#' @import ncvreg
#' @import Matrix
#' @return
#' @export
#' @param ym n*p time series
#' @param r number of factors
#' @param group  a p*r matrix, to indicate the group structure for each loading
#' @param lamvec a vector of lambda value



SparseGroupFM_bic1 <- function(ym, r, group, lamvec)
{
  p <- ncol(ym)
  n <- nrow(ym)

  Mhat <- estCov(ym, 1)
  Shat <- as.matrix(eigen(Mhat)$vectors[,1:r])

  nlam1 <- length(lamvec)
  bic1 <- rep(0, nlam1)
  loss1 <- rep(0, nlam1)
  param1 <- rep(0, nlam1)

  for(i in 1:nlam1)
  {
    resi <- SparseGroupFM(Shat = Shat, r = r, group = group,
                          lambda1 = lamvec[i], lambda2 = lamvec[i])
    if(mean(round(abs(resi$Qest),5)>0) >0)
    {
      predi <- predfun_sparse(resi,ym)
      loss1[i] <- mean((ym - predi)^2)
      param1[i] <- sum(round(resi$Qest,5)!=0)
      bic1[i] <- log(loss1[i]) + param1[i]*log(n*p)/n/p
    }
    else{break}
  }
  lammin <- lamvec[which.min(bic1)]
  res1 <- SparseGroupFM(Shat = Shat, r = r, group = group,
                        lambda1 = lammin, lambda2 = lammin)

  bic <- data.frame(lambda = lamvec, npara = param1, bic = bic1)
  out <- list(fit = res1, bic = bic, lambda = lamvec, lammin = lammin)
  return(out)
}


#'  Sparse Group Factor models based on BIC, two lambda, tuning separately. First fix lam2 as 0, and choose lam1 first, then choose lam2
#' @import MASS
#' @import ncvreg
#' @import Matrix
#' @return
#' @export
#' @param ym n*p time series
#' @param r number of factors
#' @param group  a p*r matrix, to indicate the group structure for each loading
#' @param lamvec1 a vector of lambda value for sparse
#' @param lamvec2 a vector of lambda value for group


SparseGroupFM_bic2 <- function(ym, r, group, lamvec1, lamvec2)
{
  p <- ncol(ym)
  n <- nrow(ym)

  Mhat <- estCov(ym, 1)
  Shat <- as.matrix(eigen(Mhat)$vectors[,1:r])


  ### step 1
  nlam1 <- length(lamvec1)
  bic1 <- rep(0, nlam1)
  loss1 <- rep(0, nlam1)
  param1 <- rep(0, nlam1)

  for(i in 1:nlam1)
  {
    resi <- SparseGroupFM(Shat = Shat, r = r, group = group,
                          lambda1 = lamvec1[i], lambda2 = 0)
    if(mean(round(abs(resi$Qest),5)>0) >0)
    {
      predi <- fittedfun(resi$Qest,ym)
      loss1[i] <- mean((ym - predi)^2)
      param1[i] <- sum(round(resi$Qest,5)!=0)
      bic1[i] <- log(loss1[i]) + param1[i]*log(n*p)/n/p
    }
    else{break}
  }
  lammin1 <- lamvec1[which.min(bic1)]

  res1 <- SparseGroupFM(Shat = Shat, r = r, group = group,
                        lambda1 = lammin1, lambda2 = 0)


  ### step 2

  nlam2 <- length(lamvec2)
  bic2 <- rep(0, nlam2)
  loss2 <- rep(0, nlam2)
  param2 <- rep(0, nlam2)

  for(i in 1:nlam2)
  {
    resi <- SparseGroupFM(Shat = Shat, r = r, group = group,
                          lambda1 = lammin1, lambda2 = lamvec2[i])
    if(mean(round(abs(resi$Qest),5)>0) >0)
    {
      predi <- fittedfun(resi$Qest,ym)
      loss2[i] <- mean((ym - predi)^2)
      param2[i] <- sum(round(resi$Qest,5)!=0)
      bic2[i] <- log(loss2[i]) + param2[i]*log(n*p)/n/p
    }
    else{break}
  }
  lammin2 <- lamvec2[which.min(bic2)]


  res2 <- SparseGroupFM(Shat = Shat, r = r, group = group,
                        lambda1 = lammin1, lambda2 = lammin2)

  bic1 <- data.frame(lambda = lamvec1, npara = param1, bic = bic1)
  bic2 <- data.frame(lambda = lamvec2, npara = param2, bic = bic2)
  out <- list(fit = res2,
              bic1 = bic1, bic2 = bic2,
              lambda1 = lamvec1, lambda2 = lamvec2,
              lammin1 = lammin1, lammin2 = lammin2)
  return(out)
}







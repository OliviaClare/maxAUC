#' @useDynLib maxAUC, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' @export auccp
#' ## cccp algorithm for triangular kernel using quasi-newton
# library(Rcpp)
# library(ROCR)
# library(sandwich)
# library(Metrics)
# sourceCpp('./base_func_tri_anchor_ineq.cpp')

#
#
# varauc_r <- function(coef, X, Y){
#   # variance for u statistic W
#   n = length(Y)
#   nd = sum(Y)
#   nh = n - nd
#
#   hatauc = 0
#   for(i in 1:(n-1)){
#     for(j in (i+1):n){
#       hatauc = hatauc + ifelse(X[i,]%*%coef > X[j,]%*%coef && Y[i]>Y[j], 1,0) +
#         ifelse(X[j,]%*%coef > X[i,]%*%coef && Y[j]>Y[i], 1,0)
#     }
#   }
#   # hatauc = hatauc/(n*(n-1)/2)
#   hatauc = hatauc/(nd*nh) # plug in for AUC(beta0)
#
#
#   out = 0
#   for(i in 1:n){
#     tmp = 0
#     for(j in 1:n){
#       tmp = tmp + (Y[i]>Y[j])*((X[i,]%*%coef > X[j,]%*%coef)-hatauc) +
#         (Y[j]>Y[i])*((X[j,]%*%coef > X[i,]%*%coef)-hatauc)
#     }
#     out = out + (tmp/n)^2
#   }
#   out = out/(n*(n-1))*(n*(n-1)/(nd*nh))^2
#
#   return(list("var"=out, "hatauc"=hatauc))
# }


auccp <- function(X, Y, sigma1=0.5, iter=20, outer=60, eps=1e-5, sig.eps=1e-2,
                         beta.0=NA, tol=1e-2,
                         alpha=0.8, alphat=0.8, w=2, t=1/sigma1, silence=TRUE){
  # sigma1: initial value of bandwidth
  # iter: maximum iteration for vaiance estimation
  # outer: maximum iteration for outer loop in point estimation
  # eps:

  logit <- function(x) log(x/(1-x))
  expit <- function(x) 1/(1+exp(-x))

  w.gen <- function(beta, x){
    return(exp(c(1,x)%*%beta)/(1+exp(c(1,x)%*%beta))^2)
  }

  lrvar <- function(X,Y){ #beta and var(theta) estimate from logistic regression (LR)
    d = ncol(X)-1
    b <- glm(Y~X, family = binomial())
    W = diag(apply(X, 1, w.gen, beta=b$coefficients))
    b.var = solve(t(X)%*%W%*%X)
    b.coef = b$coefficients[-1]/sum(b$coefficients[-1])

    return(list(Var=b.var, beta_initial=b$coefficients[-1]))
  }


  rbvar <- function(X,Y){
    # rb <- robustbase::BYlogreg(X,Y, addIntercept = F)
    b <- glm(Y~X, family = binomial())
    return(list(Var=sandwich(b)[-1,-1], beta_initial= b$coefficients[-1]))
  }

  delta.tf <- function(X, Y, anchor=NA, method="LR"){
    #delta transformation for LR working covariance estimate
    # anchor is from beta.est of TRI
    d = ncol(X)-1
    if(method=="LR"){
      var.est = lrvar(X,Y)$Var
      beta.est = lrvar(X,Y)$beta_initial
    }
    if(method=="RB"){
      # rb <- robustbase::BYlogreg(X,Y, addIntercept = F)
      # var.est = rb$cov
      # beta.est = rb$coefficients
      b <- glm(Y~X, family = binomial())
      var.est = sandwich(b)[-1,-1]
      beta.est = b$coefficients[-1]

    }
    # b <- glm(Y~X, family = binomial())
    # W = diag(apply(X, 1, w.gen, beta=b$coefficients))
    if(is.na(anchor)) anchor = which.max(abs(beta.est))
    # lr.var = solve(t(X)%*%W%*%X)
    # lr.beta = b$coefficients[-1]
    theta1 <- beta.est[-anchor]
    theta2 <- beta.est[anchor]
    sig11 = var.est[-anchor, -anchor]
    sig22 = var.est[anchor, anchor]
    sig12 = var.est[-anchor, anchor]
    out =sig11/theta2^2 - 2*theta1%*%t(sig12)/theta2^3 + sig22*theta1%*%t(theta1)/theta2^4
    return(out)
  }

  warn_msgs <- character()
  suppressWarnings({
    # Your code that generates warnings

  # withCallingHandlers(
  #   {



  # triangular kernel based CCP optimization
  X <- as.matrix(X)
  Y <- as.vector(Y)
  d = ncol(X)-1
  n = length(Y)
  Xd = X[which(Y==1),]
  Xh = X[which(Y==0),]
  nd = nrow(Xd)
  nh = nrow(Xh)
  Y <- c(rep(1,nd),rep(0,nh))
  X = rbind(Xd,Xh)



  # the objective function to be optimized
  find.lam.kernel <- function(beta){
    -tauc_sort(beta, X, Y, sigma1, T) - w*sum(beta^2) - log(1-sum(beta^2))/t
  }


  # the gradient function
  grn <- function(beta){
    dtauc_opt(beta, beta.k, X, Y, sigma1, 2, t)
  }




  #step 0
  if(is.na(sum(beta.0))){
    b <- glm(Y~X, family = binomial())
    if(sum(b$coefficients[-1]<0)>0){
      mu_delta.sign <- sign(b$coefficients[-1])
      X = t(apply(X,1,function(x) x*mu_delta.sign))
      b <- glm(Y~X, family = binomial())
    }else{mu_delta.sign <- sign(b$coefficients[-1])}

    beta.k <- b$coefficients[-1]/norm(b$coefficients[-1], type="2")
  }
  else {
    beta.k = beta.0
    mu_delta.sign <- sign(beta.k)
    X = t(apply(X,1,function(x) x* mu_delta.sign))
  }

  beta.k = beta.k*0.8

  p = ncol(Xd)
  ptm <- proc.time()
  outer.history = matrix(rep(NA,(1+p)*outer),nrow=outer)
  # step 1 to k
  for(k in 1:outer){
    beta.in = optim(beta.k, find.lam.kernel, method="BFGS", gr=grn,
                          control=list(maxit=30000, ndeps=.1))$par
    eauc.current = eauc_sort(beta.in, X, Y, silence = TRUE)
    # cat("l2norm", sqrt(sum(beta.in^2)),'\n')
    # beta.in = beta.in*0.8

    outer.history[k,] = c(eauc.current, beta.in)
    if(silence==FALSE) cat(outer.history[k,],'\n')
    # if(eauc_sort(beta.1,X,Y,silence=TRUE)-eauc_sort(beta.0,X,Y,silence=TRUE)<eps) beta.1 = 0.5*(beta.1+beta.0)
    sigma1 = sigma1*alpha
    t = t/alphat
    beta.k = beta.in
    if(sigma1 < sig.eps){
      if(abs(outer.history[k,1]-outer.history[k-1,1])<eps) break
    }
  }
  if(!silence) cat("eauc history:", outer.history[1:k,1])

  #   },
  # warning = function(w) warn_msgs <<- append(warn_msgs, w$message)
  # )
  })
  result = outer.history[which(outer.history[,1]==max(outer.history[,1],na.rm = T))[1], ]
  beta.est = result[1+1:p]/sqrt(sum(result[1+1:p]^2))
  time.use = sum((proc.time()-ptm)[1:2])

  # inference based on Zhang et al.
  anchor = which.max(abs(beta.est))
  # theta.est = beta.est/abs(beta.est[anchor])
  cov.es1 = dn(beta.est,X,Y,delta.tf(X,Y,anchor)*n, anchor)

  ###############################################################################
  # -----------------------------------------------------------------------------
  for(i in 1:iter){
    # z.var <- dn_anchor(theta.est,X,Y,as.matrix(z.var))*n
    cov.es.new = dn(beta.est,X,Y,cov.es1, anchor)
    if(is.na(sum(cov.es.new))){
      if(!silence) cat('did not converge.\n')
      cov.es1 = delta.tf(X,Y,anchor)*n ## times n
      break
    }
    if(!silence) cat("new cov:", cov.es.new, "\n")
    if(sum(abs(cov.es.new-cov.es1)) < tol) {if(!silence) cat('converge to:',cov.es.new,'\n'); break}
    cov.es1 = cov.es.new
  }

  #------------------------------------------------------------------------------------
  ##############################################################################

  hessian.es1 = an(beta.est, X,Y, cov.es1, anchor)
  meat = vn(beta.est, X,Y,cov.es1, anchor)


  # bias correction
  Bias = -sum(diag(hessian.es1%*%cov.es1/n))
  bias_corrected_pAUC = result[1] - Bias # bias corrected predictive AUC (AUC(\hatbeta))
  bias_corrected_tAUC = result[1] - Bias/2 # bias corrected theoretical AUC (AUC(beta0))

  # logit bias correction
  bias.logit = (1/result[1]+1/(1-result[1]))*Bias
  bias_corrected_logit_pAUC = expit(logit(result[1])-bias.logit)
  bias_corrected_logit_tAUC = expit(logit(result[1])-bias.logit/2)


  # variance of AUC
  var.auc= varauc_l1(beta.est, X,Y)
  # var.auc = (n-1)^2/(4*n^2)*(nd*nh)^2/n^4*varauc_l1(beta.est,X,Y) + (n-1)^2/(4*n^3)*(1-2*pi.est)^2*pi.est*(1-pi.est)*result[1]
  CIlogit <- function(point.est, var.auc){
    logit.est = logit(point.est)
    varlogitAUC = var.auc*(1/point.est+1/(1-point.est))^2
    left = expit(logit.est - 1.96*sqrt(varlogitAUC))
    right = expit(logit.est + 1.96*sqrt(varlogitAUC))
    return(c(point.est, left, right))
  }
  bias_correction <- rbind(c(bias_corrected_pAUC,
                             bias_corrected_pAUC-1.96*sqrt(var.auc),
                             bias_corrected_pAUC+1.96*sqrt(var.auc)),
                           c(bias_corrected_tAUC,
                             bias_corrected_tAUC-1.96*sqrt(var.auc),
                             bias_corrected_tAUC+1.96*sqrt(var.auc)),
                           CIlogit(bias_corrected_logit_pAUC, var.auc),
                           CIlogit(bias_corrected_logit_tAUC, var.auc))
  bias_correction <- as.data.frame(bias_correction)
  colnames(bias_correction) = c("AUC", "CI.lower","CI.upper")
  rownames(bias_correction) = c("bias_corrected_pAUC", "bias_corrected_tAUC",
                                "bias_corrected_logit_pAUC", "bias_corrected_logit_tAUC")

  return(list("eAUC"=result[1], "beta"=mu_delta.sign*beta.est,"time.elapsed"= time.use, "var.auc"=var.auc,
              "bias_correction"=bias_correction,
              "cov"=cov.es1))
  # , nwarnings = length(warn_msgs)
}

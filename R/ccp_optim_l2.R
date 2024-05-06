#' @useDynLib maxAUC, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' @export aucmax
#' ## cccp algorithm for triangular kernel using quasi-newton



## cccp algorithm for combination estimation with variance estimation procedure

# library(Rcpp)
# library(Metrics)


# # import functions: eauc_sort; tauc_sort; dtauc_opt; varauc
# # import functions: varbeta; hessianpw; meatpw

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))

lrvar <- function(X,Y){ #beta and var(theta) estimate from logistic regression (LR)

  # w.gen <- function(beta, x){
  #   return(exp(c(1,x)%*%beta)/(1+exp(c(1,x)%*%beta))^2)
  # }

  b <- glm(Y~X, family = binomial())
  # W = diag(apply(X, 1, w.gen, beta=b$coefficients))
  # b.var = solve(t(X)%*%W%*%X)
  didv = norm(b$coefficients[-1],"2")
  b.coef = b$coefficients[-1]/didv
  b.cov = vcov(b)[-1,-1]/didv^2
  return(list(Var=b.cov, beta_initial=b$coefficients[-1]/didv))
}

aucmax <- function(X, Y, sigma1=1, sig.eps=1e-2, outer=60, eps=1e-8,
                            beta.0=NA, alpha=0.9, w=2,
                            var.method="selfinduced", tol=2, inf_iter = 10,
                            silence=FALSE){
  # sigma1: initial value of bandwidth
  # sig.eps: upper bound for sigma1 for convergence
  # eps: upper bound of convergence for empirical AUC in combination estimation
  # outer: maximum iteration for outer loop in point estimation
  # beta.0: initial beta estimate
  # alpha: shrinkage factor for bandwidth
  # w: penalty parameter for l2 constraint

  # var.method: variance estimation method, can be numdrv or selfinduced
  # tol: upper bound of convergence for asymptotic variance
  # inf_iter: maximum iteration for variance estimation


  warn_msgs <- character()
  suppressWarnings({
    # code that generates warnings

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
    wobj <- function(beta){
      -tauc_sort(beta, X, Y, sigma1, T) +
        w*(1-norm(beta,"2"))^2
      # w*(sum(beta^2) - 2*sum(beta.k*beta)/sqrt(sum(beta.k^2)))
    }


    # the gradient function
    grn <- function(beta){
      dtauc_opt(beta, beta.k, X, Y, sigma1) + 2*w *
        (beta-beta.k/sqrt(sum(beta.k^2)))
    }

    #step 0
    if(is.na(sum(beta.0))){
      b <- glm(Y~X, family = binomial())
      if(is.na(sum(b$coefficients))) {
        beta.k = rep(1, d+1)/norm(rep(1, d+1),"2")
        mu_delta.sign = 1
      } else{
        if(sum(b$coefficients[-1]<0)>0){
          mu_delta.sign <- sign(b$coefficients[-1])
          X = t(apply(X,1,function(x) x*mu_delta.sign))
          # b <- glm(Y~X, family = binomial())
        }else{mu_delta.sign <- sign(b$coefficients[-1])}

        beta.k <- b$coefficients[-1]/norm(b$coefficients[-1], type="2")
      }
    }
    else {
      mu_delta.sign <- sign(beta.0)
      X = t(apply(X,1,function(x) x* mu_delta.sign))
      beta.k = abs(beta.0)
    }

    # beta.k = beta.k*0.8

    p = ncol(Xd)
    ptm <- proc.time()
    outer.history = matrix(rep(NA,(1+p)*outer),nrow=outer)
    # step 1 to k
    for(k in 1:outer){
      eauc.current = eauc_sort(beta.k, X, Y, silence = TRUE)
      # cat("l2norm", sqrt(sum(beta.in^2)),'\n')
      # beta.in = beta.in*0.8

      outer.history[k,] = c(eauc.current, beta.k)

      beta.in = optim(beta.k, wobj, method="BFGS", gr=grn,
                      control=list(maxit=30000, ndeps=.1))$par

      # if(silence==FALSE) cat(outer.history[k,],"w=",w, '\n')
      # if(eauc_sort(beta.1,X,Y,silence=TRUE)-eauc_sort(beta.0,X,Y,silence=TRUE)<eps) beta.1 = 0.5*(beta.1+beta.0)
      sigma1 = sigma1*alpha
      # t = t/alphat
      # w = w/alphaw
      beta.k = beta.in
      if(sigma1 < sig.eps){
        if(abs(outer.history[k,1]-outer.history[k-1,1])<eps) break
      }
    }
    # cat("eauc history:", outer.history[1:k,1])

    #   },
    # warning = function(w) warn_msgs <<- append(warn_msgs, w$message)
    # )
  })
  result = outer.history[which(outer.history[,1]==max(outer.history[,1],na.rm = T))[1], ]
  beta.est = result[1+1:p]/sqrt(sum(result[1+1:p]^2))
  time.use = sum((proc.time()-ptm)[1:2])

  ###############################################################################
  # variance estimation -----------------------------------------------------------------------------
  if(is.na(var.method)){
    # cov.es1 = NA
    return(list("eAUC"=result[1], "beta"=mu_delta.sign*beta.est,
                "time.elapsed"= time.use,
                # "var.auc"=var.auc,
                # "bias_correction"=bias_correction,
                "cov"=NA
                # "cov.iter"=cov.re$iter,
                # "cov.converge" = cov.re$converge
    ))}
  else if(var.method=="selfinduced"){
    var.ptm <- proc.time()
    # inference based on Zhang et al.
    # anchor = which.max(abs(beta.est))
    # theta.est = beta.est/abs(beta.est[anchor])
    cov.es1 = varbeta(beta.est,X,Y,lrvar(X,Y)$Var*n, w)

    for(i in 1:inf_iter){
      # z.var <- dn_anchor(theta.est,X,Y,as.matrix(z.var))*n
      cov.es.new = varbeta(beta.est,X,Y,cov.es1, w)
      if(is.na(sum(cov.es.new))){
        if(!silence) cat('did not converge.\n')
        cov.es1 = lrvar(X,Y)*n ## times n
        break
      }
      if(!silence) {
        cat("new meat:", meatpw(beta.est, X,Y,cov.es1)*n^3/(nd^2*nh^2),"\n")
        cat("new cov:", cov.es.new, "\n")}
      if(sum(abs(cov.es.new-cov.es1)) < tol) {if(!silence) cat('converge to:',cov.es.new,'\n'); break}
      cov.es1 = cov.es.new
    }

    hessian.es1 = hessianpw(beta.est, X,Y, cov.es1, w)*n*(n-1)/(nd*nh) #adjust to auc scale
    # meat = meatpw(beta.est, X,Y,cov.es1)

    var.time.use = sum((proc.time()-var.ptm)[1:2])

    # bias correction
    Bias = -sum(diag(hessian.es1%*%cov.es1/n))
    bias_corrected_pAUC = result[1] - Bias # bias corrected predictive AUC (AUC(\hatbeta))
    bias_corrected_tAUC = result[1] - Bias/2 # bias corrected theoretical AUC (AUC(beta0))
    # AUC_bic = result[1]- Bias*log(n)/2

    # logit bias correction
    bias.logit = (1/result[1]+1/(1-result[1]))*Bias
    bias_corrected_logit_pAUC = expit(logit(result[1])-bias.logit)
    bias_corrected_logit_tAUC = expit(logit(result[1])-bias.logit/2)


    # variance of AUC
    var.auc= varauc(beta.est, X,Y)
    # var.auc = (n-1)^2/(4*n^2)*(nd*nh)^2/n^4*varauc(beta.est,X,Y) + (n-1)^2/(4*n^3)*(1-2*pi.est)^2*pi.est*(1-pi.est)*result[1]
    #
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

    # , nwarnings = length(warn_msgs)

    return(list("eAUC"=result[1], "beta"=mu_delta.sign*beta.est,
                "time.elapsed"= time.use,
                "var.time"=var.time.use,
                "var.auc"=var.auc,
                "bias_correction"=bias_correction,
                # "bic"=AUC_bic,
                "cov"=cov.es1/n
    ))
  }

}

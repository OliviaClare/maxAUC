#' @useDynLib maxAUC, .registration = TRUE
#' @importFrom Rcpp sourceCpp

#' iterative algorithm to find theta estimate
#' @export smrce

smrce <- function(X,Y, beta_initial, tol,iteration, silence = TRUE){

  w.gen <- function(beta, x){
    return(exp(c(1,x)%*%beta)/(1+exp(c(1,x)%*%beta))^2)
  }

  lrvar <- function(X,Y){ #beta and var(theta) estimate from logistic regression (LR)
    d = ncol(X)-1
    b <- glm(Y~X, family = binomial())
    W = diag(apply(X, 1, w.gen, beta=b$coefficients))
    b.var = solve(t(X)%*%W%*%X)
    b.coef = b$coefficients[-1]/sum(b$coefficients[-1])

    return(list(Var=b.var, beta_initial=b.coef))
  }


  converge=1
  ptm <- proc.time()
  # make sure the input format is correct
  X <- as.matrix(X)
  Y <- as.vector(Y)
  d = ncol(X)-1
  n = length(Y)
  Xd = X[which(Y==1),]
  Xh = X[which(Y==0),]
  X_delta = colMeans(Xd)-colMeans(Xh)

  auc.score =  pnorm(X_delta/sqrt(2))
  bench = d+1 #default bench biomarker (single biomarker with max AUC)
  if(auc.score[d+1] != max(auc.score)){ #set bench biomarker
    if(!silence) cat("Warning: INCORRECT ANCHOR BIOMARKER!\n")
    bench = which(auc.score==max(auc.score))
    tmp = Xd[,d+1]
    Xd[,d+1] = Xd[,bench]
    Xd[,bench] = tmp
    tmp = Xh[,d+1]
    Xh[,d+1] = Xh[,bench]
    Xh[,bench] = tmp
    theta_initial = beta_initial[1:d]
    theta_initial[bench] = beta_initial[d+1]
    theta_initial = theta_initial/abs(beta_initial[bench])
  }
  else{
    theta_initial = beta_initial[1:d]/abs(beta_initial[bench])
  }
  X = rbind(Xd,Xh)
  Y = c(rep(1,nrow(Xd)),rep(0,nrow(Xh)))


  #initialize variance of z
  z.var <- n*lrvar(X,Y)$Var[1:d,1:d]
  theta.max = theta.new = theta_initial
  qn.max = qn.old = eauc(theta.max,X,Y)

  for (k in 1:iteration){
    # if(sum(abs(theta.max))>4*d) { #reinitiate if the estimate is too deviated
    #   theta.max = runif(d,0,1)
    #   z.var <- 100*diag(d)
    # }
    z.var <- dn_anchor(theta.new,X,Y,as.matrix(z.var)) #new covariate matrix
    tmp <- try(newton_raphson_anchor(theta.max,tol=tol, iteration = 50, gamma=0.5, X,Y,z.var), silent = T)
    if(!is.numeric(tmp)) {
      converge=0
      break
    }
    else theta.new = tmp
    qn.new = eauc(theta.new,X,Y) #corresponding value of objective function

    if(qn.new > qn.max) {qn.max = qn.new;     theta.max = theta.new}
    if(!silence) cat('iteration:',k,'\t theta:',theta.new,'\t', 'qn:', qn.new,'\n', 'var:',z.var,'\n')
    if(abs(qn.old-qn.new) < tol && k>4) {if(!silence) cat('converge to:',theta.max,'\n'); break}
    qn.old = qn.new
  }

  # ###############if not converging##################
  # if(k == iteration){
  #   cat("Warning: Algorithm did not converge.")
  #   #brutal grid search
  #   theta.grid=do.call(expand.grid, replicate(d, seq(0,2,0.1), simplify=FALSE))
  #   z.var=20*diag(d)
  #   for (i in 1:nrow(theta.grid)){
  #     theta.grid[i,d+1] = eauc(as.double(theta.grid[i,1:d]),X,Y)
  #   }
  #   theta.max = as.numeric(theta.grid[which(theta.grid[,d+1]==max(theta.grid[,d+1]))[1],1:d])
  #   z.var = matrix(0,d,d)
  # }

  eAUC = eauc(theta.max,X,Y)
  # sAUC = ifelse(k==iteration, 0, qncpp(theta.max,X,Y,z.var))
  hessian.es = an_anchor(theta.max,X,Y,z.var)
  meat.es = vn_anchor(theta.max,X,Y,z.var)
  if(bench != d+1){ #reparameterize if bench biomarker is not default
    theta.oldbench = theta.max[bench]
    theta.max[bench] = 1
    beta.max = c(theta.max,theta.oldbench)
  }
  else beta.max = c(theta.max, 1)
  #use l2 norm to constrain beta
  beta.max = beta.max/norm(beta.max, type="2")


  time.elapsed = sum((proc.time()-ptm)[1:2])
  return_list = list("eAUC"=eAUC, "beta"=beta.max, "cov"=z.var,
                     "iter"=k,
                     "time.elapsed"=time.elapsed, "converged"= converge)
  return(return_list)
}


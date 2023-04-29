#' @useDynLib maxAUC, .registration = TRUE

#' @export sim4
sim4 <- function(sample.size=100, scenario=1, markers=3, pd=0.3, design="cohort"){
  #simulate dataset
  # ph = 0.7 #number of controls
  # pd = 0.3 #number of cases
  ph = 1 - pd #number of controls

  if(scenario==1){ # 0.8648279
    mud = rep(0.9, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    cov.tr1 = diag(rep(1,markers))

  }else if(scenario==2){ # 0.8364066
    mud = rep(0.8, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    if(markers==3){
      cov.tr1 = diag(c(0.5, 1, 2))
    }else if(markers==6){
      cov.tr1 = diag(c(0.5, 1, 2, 1, 1, 1))
    }
  }
  else if(scenario==3){ # 0.843312
    mud = rep(1, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    if(markers==3){
      cov.tr1 = diag(c(0.5, 1, 2))
      cov.tr1[which(cov.tr1==0)]=0.5
    }else if(markers==6){
      cov.tr1 = diag(c(0.5, 1, 2, 1, 1, 1))
      cov.tr1[1:3,1:3] = 0.5
      cov.tr1[2,2] = 1
      cov.tr1[3,3] = 2
    }
  }
  else if(scenario==4){  # AUC.tr = 0.8636912
    p.subtype = 2/3 # prob of subtypes
    mud1 = c(1.7, 1.7, 0) #mean of subtype1 group
    mud2 = c(0, 0, 1.7) #mean of subtype2 group
    # mud3 = c(0.5,0.5) #mean of subtype3 group
    muh = rep(0, markers) #mean of control group
    cov.trd1 = diag(c(0.5,2,1))
    cov.trd2 = diag(rep(1,markers))
    if(markers==6){
      mud1 = c(mud1, 0, 0, 0)
      mud2 = c(mud2, 0, 0, 0)
      cov.trd1 = diag(c(0.5,2,1,1,1,1))
    }
  }

  d = markers-1
  cov.tr0 = diag(rep(1, markers))
  # mu_delta = mud-muh
  # beta.true = mu_delta%*%solve(cov.tr1)/sum(mu_delta%*%solve(cov.tr1)) #theoretical coefficient
  # beta.true = beta.true/sum(abs(beta.true)) #contrain by l2 norm
  # theta.true <- beta.true[1:d]/sum(beta.true)
  # #theoratical AUC
  # AUC.tr = pnorm(sum(beta.true*mu_delta)/sqrt(beta.true%*%(cov.tr1+cov.tr0)%*%t(beta.true)))
  # AUC.tr = max(AUC.tr, 1-AUC.tr)


  #simulate training data
  if(scenario %in% 1:3){
    if(design=="cohort"){
      nd = rbinom(1, sample.size, pd)
      nh = sample.size-nd
    }else if(design=="casecontrol"){
      nd = sample.size*pd
      nh = sample.size*ph
    }
    xd = mvrnorm(nd,mud, cov.tr1)
    xh = mvrnorm(nh,muh, cov.tr0)
  }
  if(scenario == 4){
    #simulate training data
    if(design=="cohort"){
      nd = rbinom(1,sample.size,pd)
      nd1 = rbinom(1,nd,p.subtype)
      nd2 = nd - nd1
      nh = sample.size - sum(nd)
    }else if(design=="casecontrol"){
      nd1 = sample.size*pd*p.subtype
      nd2 = sample.size*pd - nd1
      nh = sample.size*ph
    }
    xd = rbind(mvrnorm(nd1,mud1, cov.trd1),mvrnorm(nd2,mud2, cov.trd2))
    xh = mvrnorm(nh,muh, cov.tr0)
  }

  # xh[c(1:20),] = rep(10, (d+1))
  # x_delta = colMeans(xd)-colMeans(xh)
  Y <- c(rep(1,nd),rep(0,nh))
  X = rbind(xd,xh)
  colnames(X) = paste("x",1:markers, sep="")
  # data set
  res <- data.frame(y=Y, X)


}

tAUC <- function(beta, scenario=4){
  markers = length(beta)
  if(scenario==1){ # 0.8648279
    mud = rep(0.9, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    cov.tr1 = diag(rep(1,markers))

  }else if(scenario==2){ # 0.855040773
    mud = rep(0.8, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    if(markers==3){
      cov.tr1 = diag(c(0.5, 1, 2))
    }else if(markers==6){
      cov.tr1 = diag(c(0.5, 1, 2, 1, 1, 1))
    }
  }
  else if(scenario==3){ # 0.841344746
    mud = rep(1, markers)
    if(markers==6) mud[4:6] = 0
    muh = rep(0, markers)
    rho=0
    if(markers==3){
      cov.tr1 = diag(c(0.5, 1, 2))
      cov.tr1[which(cov.tr1==0)]=0.5
    }else if(markers==6){
      cov.tr1 = diag(c(0.5, 1, 2, 1, 1, 1))
      cov.tr1[1:3,1:3] = 0.5
      cov.tr1[2,2] = 1
      cov.tr1[3,3] = 2
    }
  }else if(scenario==4){  # AUC.tr = 0.8636912
    p.subtype = 2/3 # prob of subtypes
    mud1 = c(1.7, 1.7, 0) #mean of subtype1 group
    mud2 = c(0, 0, 1.7) #mean of subtype2 group
    # mud3 = c(0.5,0.5) #mean of subtype3 group
    muh = rep(0, markers) #mean of control group
    cov.trd1 = diag(c(0.5,2,1))
    cov.trd2 = diag(rep(1,markers))
    if(markers==6){
      mud1 = c(mud1, 0, 0, 0)
      mud2 = c(mud2, 0, 0, 0)
      cov.trd1 = diag(c(0.5,2,1,1,1,1))
    }
  }

  d = markers-1
  cov.tr0 = diag(rep(1, markers))

  if(scenario %in% 1:3){
      mu_delta = mud-muh
    #theoratical AUC
      AUC = pnorm(sum(beta*mu_delta)/sqrt(beta%*%(cov.tr1+cov.tr0)%*%beta))
  }
  if(scenario == 4){
    AUC=p.subtype*pnorm(sum(beta*mud1)/sqrt(beta%*%(cov.trd1+cov.tr0)%*%beta)) +
      (1-p.subtype)*pnorm(sum(beta*mud2)/sqrt(beta%*%(cov.trd2+cov.tr0)%*%beta))

  }
  AUC = max(AUC, 1-AUC)
  return(as.numeric(AUC))
}


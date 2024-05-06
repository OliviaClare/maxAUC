#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

// [[Rcpp::export]]
arma::vec ddnormcpp(arma::vec x){
  arma::vec out;
  out = -(x*arma::normpdf(x));
  return out;
}



// [[Rcpp::export]]
arma::mat hessianpw(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, double w){
  // calculates hessian
  // beta should have anchor set to 1
  // beta = beta/beta[anchor-1];
  int d = X.n_cols;
  arma::mat out; out=arma::zeros(d,d);
  arma::vec tmp;
  // arma::vec beta = add_one(theta,1,anchor);
  // arma::uvec v = find(beta!=1); //index vector for non-anchor biomarkers
  // v.print("v=");
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;
  int n = nd + nh;
  arma::mat mtmp; mtmp=arma::zeros(d,d);

  for(int i=0; i<nd; i++){
    for(int j=0; j<nh; j++){
      arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]); //row vector 1*(d+1)
      // tmp_x.print("tmp_x");
      // tmp_x.elem(v).print("tmp_x.v");
      arma::mat tmp_sig = sqrt(tmp_x*var*tmp_x.t());
      // tmp_sig.print("tmp_sig");
      if(tmp_sig.is_zero()){
        continue;
      }
      mtmp = tmp_x.t()*(n*ddnormcpp(sqrt(n)*tmp_x*beta/tmp_sig)/square(tmp_sig))*tmp_x;
      out = out + mtmp;
    }

  }

  out = out/(n*(n-1));
  // out = out/(nd*nh);
  //correct Hessian to semi-negative definite
  arma::vec eigval;
  arma::mat eigvec;
  eig_sym(eigval, eigvec, out);
  // eigval.print("eigval");
  // eigvec.print("eigvec");
  // out.print("out");
  out = out - w*beta*beta.t()/pow(norm(beta),3) - w*(norm(beta)-1)/norm(beta) * arma::eye<arma::fmat>(d,d);
  // eig_sym(eigval, eigvec, out);
  // eigval.elem(find(eigval>0)).zeros();
  // eigval.elem(find(abs(eigval)<1e-4)).zeros();
  // eigval.print("eigval");
  // eigvec.print("eigvec");
  // // eigval.print("eigval");
  // // out = eigvec*diagmat(eigval)*eigvec.t();
  // out.print("out2");
  // if((-out).is_sympd()){
  // 	// printf("neg def \n");
  // }else{
  // 	// out.print("out=");

  // }


  return out;
}

// [[Rcpp::export]]
arma::mat meatpw(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var){
  // beta should have anchor set to 1
  // beta = beta/beta[anchor-1];
  int d = X.n_cols;
  arma::mat out; out=arma::zeros(d,d);
  arma::vec tmp;
  // arma::vec beta = add_one(theta,1,anchor);
  // arma::uvec v = find(beta!=1); //index vector for non-anchor biomarkers
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;
  int n = nd + nh;

  for(int i=0; i<nd; i++){
    tmp.zeros(d);
    for(int j=0; j<nh; j++){
      arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]); //row vector 1*(d+1)
      arma::mat tmp_sig = sqrt(tmp_x*var*tmp_x.t());
      if(tmp_x.is_zero()){
        // out += 0.5;
        continue;
      }
      // tmp.print("tmp=");
      tmp = tmp + tmp_x.t()*(normpdf(sqrt(n)*tmp_x*beta/tmp_sig)*sqrt(n)/tmp_sig);
    }
    out = out + tmp*tmp.t();
    // out.print("out=");
  }

  for(int i=0; i<nh; i++){
    tmp.zeros();
    for(int j=0; j<nd; j++){
      arma::mat tmp_x = X.row(hinx[i])-X.row(dinx[j]); //row vector 1*(d+1)
      arma::mat tmp_sig = sqrt(tmp_x*var*tmp_x.t());
      if(tmp_x.is_zero()){
        // out += 0.5;
        continue;
      }
      tmp = tmp + tmp_x.t()*(normpdf(sqrt(n)*tmp_x*beta/tmp_sig)*sqrt(n)/tmp_sig);
    }
    out = out + tmp*tmp.t();
  }

  out = out/pow(n,3);
  // out = out*n/(pow(nd,2)*pow(nh,2));
  return out;
}



// [[Rcpp::export]]
arma::mat varbeta(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var,double w){
  // int n = X.n_rows;
  arma::mat an_inv;
  arma::mat out;
  // beta = beta/beta[anchor-1];
  // beta.print("beta=");
  an_inv = pinv(hessianpw(beta,X,Y,var,w));
  out = an_inv*meatpw(beta,X,Y,var)*an_inv;
  return out;
}

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

// [[Rcpp::export]]
double eauc_sort(arma::vec beta, arma::mat X, arma::vec Y, bool silence){
  double out = 0;
  int d = X.n_cols - 1;
  // arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
  // arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
  int beta_d = beta.n_elem;
  if(beta_d!=(d+1)){
    printf("The dimension of coefficient is incorrect.\n");
    return 0;
  }
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;
  //hinx.print("hinx=");

  arma::mat xd = X.rows(dinx);
  arma::mat xh = X.rows(hinx);

  // calculate and sort the linear combination for each group
  arma::vec dlc = reverse(sort(xd*beta)); //default is ascending: sort_direction="ascend"
  arma::vec hlc = sort(xh*beta);
  // dlc.print("dlc=");
  // xd.print("xd");


  for (int i = 0; i < nd; i++){
    if(dlc[i]-hlc[0]<0) continue;
    for  (int j = 0; j <nh; j++){
      if(dlc[i]-hlc[j] > 0){
        out += 1;
      }else continue;
    }
  }

  out = out/(nd*nh);

  out = (out>=0.5)*out + (out<0.5)*(1-out);

  return out;
}

double triang(double x, double sigma0){
  double out;
  if(x < -sigma0){
    out = 0;
    return out;
  }else if(x < 0){
    out = x*x/(2*sigma0*sigma0) + x/sigma0 + 0.5;
    return out;
  }else if(x <sigma0){
    out = -x*x/(2*sigma0*sigma0) + x/sigma0 + 0.5;
    return out;
  }else{
    out = 1;
    return out;
  }
}

// [[Rcpp::export]]
double tauc_sort(arma::vec beta, arma::mat X, arma::vec Y, double sigma0, bool silence){
  double out = 0;
  int d = X.n_cols - 1;
  // arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
  // arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
  int beta_d = beta.n_elem;
  if(beta_d!=(d+1)){
    printf("The dimension of coefficient is incorrect.\n");
    return 0;
  }
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;

  arma::mat xd = X.rows(dinx);
  arma::mat xh = X.rows(hinx);

  // calculate and sort the linear combination for each group
  arma::vec dlc = reverse(sort(xd*beta)); //default is ascending: sort_direction="ascend"
  arma::vec hlc = sort(xh*beta);


  for (int i = 0; i < nd; i++){
    if(dlc[i]-hlc[0]<-sigma0) continue;
    for  (int j = 0; j <nh; j++){
      double tmp = dlc[i]-hlc[j];
      if(tmp > -sigma0){
        out += triang(tmp, sigma0);
      }else continue;
    }
  }


  out = out/(nd*nh);

  // out = out;

  return out;
}
//
// // [[Rcpp::export]]
// double tauc_sort_cccp(arma::vec beta, arma::vec beta_k, arma::mat X, arma::vec Y, double sigma0, bool silence){
//   double out = 0;
//   // int d = X.n_cols;
//
//   // arma::vec out; out=arma::zeros(d);
//   arma::uvec dinx = find(Y==1);
//   arma::uvec hinx = find(Y==0);
//   arma::mat delta_k;
//   arma::mat delta;
//   int nd = dinx.n_elem;
//   int nh = hinx.n_elem;
//   arma::mat xd = X.rows(dinx);
//   arma::mat xh = X.rows(hinx);
//
//   for (int i = 0; i < nd; i++){
//     // xd.row(i).print("row i ");
//     arma::mat tmp_x = -(xh.each_row() - xd.row(i)); //add to each row of xh
//     // tmp_x.print("tmp_x=");
//     arma::mat delta = sum(tmp_x.each_row()%(beta.t()),1);
//     // delta.print("delta=");
//     arma::mat delta_k = sum(tmp_x.each_row()%(beta_k.t()),1);
//     for (int j = 0; j < nh; j++){
//       if(delta_k[j] + sigma0 < 0){
//         if(delta[j] < 0) continue;
//         else if(delta[j] < sigma0) out = out + pow(delta[j],2)/(2*pow(sigma0,2))+delta[j]/sigma0;
//         else out = out + 2*delta[j]/sigma0 - 0.5;
//       }
//       else if(delta_k[j] < 0){
//         // delta.print("delta=");
//         // tmp_x.print("tmp_x=");
//         // out.print("out=");
//         if(delta[j] < 0) out = out - delta[j]/sigma0 + (delta_k[j]*(delta_k[j]-2*delta[j]))/(2*pow(sigma0,2)) - 0.5;
//         else if(delta[j] < sigma0) out = out + pow((delta[j]-delta_k[j])/sigma0,2)/2 - 0.5;
//         else out = out + delta[j]/sigma0 + (delta_k[j]*(delta_k[j]-2*delta[j]))/(2*pow(sigma0,2)) - 1;
//       }
//       else{
//         // delta.print("delta=");
//         if(delta[j] < 0) out = out - 2*delta[j]/sigma0 - 0.5;
//         else if(delta[j] < sigma0) out = out  + pow(delta[j],2)/(2*pow(sigma0,2)) - delta[j]/sigma0;
//         else out = out + 1;
//       }
//       // out.print("out=");
//     }
//   }
//   out = out/(nd*nh);
//   // out = out - w * 2 * beta_k;
//   // out.print("out="); beta_k.print("beta=");
//   // out = out + 2*beta/(t*(1-pow(norm(beta),2)));
//   return out;
// }



// [[Rcpp::export]]
arma::mat dtauc_opt(arma::vec beta, arma::vec beta_k, arma::mat X, arma::vec Y, double sigma0){
  // to be minimized
  int d = X.n_cols;

  arma::vec out; out=arma::zeros(d);
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  arma::mat delta_k;
  arma::mat delta;
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;
  arma::mat xd = X.rows(dinx);
  arma::mat xh = X.rows(hinx);

  for (int i = 0; i < nd; i++){
    // xd.row(i).print("row i ");
    arma::mat tmp_x = -(xh.each_row() - xd.row(i)); //add to each row of xh
    // tmp_x.print("tmp_x=");
    arma::mat delta = sum(tmp_x.each_row()%(beta.t()),1);
    // delta.print("delta=");
    arma::mat delta_k = sum(tmp_x.each_row()%(beta_k.t()),1);
    for (int j = 0; j < nh; j++){
      if(delta_k[j] + sigma0 < 0){
        if(delta[j] < 0) continue;
        else if(delta[j] < sigma0) out = out + (delta[j]/pow(sigma0,2)+1/sigma0)*tmp_x.row(j).t();
        else out = out + 2*tmp_x.row(j).t()/sigma0;
      }
      else if(delta_k[j] < 0){
        // delta.print("delta=");
        // tmp_x.print("tmp_x=");
        // out.print("out=");
        if(delta[j] < 0) out = out - (delta_k[j]/pow(sigma0,2)+1/sigma0)*tmp_x.row(j).t();
        else if(delta[j] < sigma0) out = out + (delta[j]-delta_k[j])/pow(sigma0,2)*tmp_x.row(j).t();
        else out = out + (1/sigma0-delta_k[j]/pow(sigma0,2))*tmp_x.row(j).t();
      }
      else{
        // delta.print("delta=");
        if(delta[j] < 0) out = out - 2*tmp_x.row(j).t()/sigma0;
        else if(delta[j] < sigma0) out = out + (delta[j]/pow(sigma0,2)-1/sigma0)*tmp_x.row(j).t();
        else continue;
      }
      // out.print("out=");
    }
  }
  out = out/(nd*nh);
  // out = out - w * 2 * beta_k;
  // out.print("out="); beta_k.print("beta=");
  // out = out + 2*beta/(t*(1-pow(norm(beta),2)));
  return out;
}




// [[Rcpp::export]]
double varauc(arma::vec beta, arma::mat X, arma::vec Y){
  double eauc = 0;
  double out1 = 0; double out2 = 0;
  int d = X.n_cols - 1;

  int beta_d = beta.n_elem;
  if(beta_d!=(d+1)){
    printf("The dimension of coefficient is incorrect.\n");
    return 0;
  }
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;

  for (int i = 0; i < nd; i++){
    for (int j = 0; j < nh; j++){
      if(arma::as_scalar((X.row(dinx[i])-X.row(hinx[j]))*beta)>0) eauc += 1;
    }
  }

  eauc = eauc/(nd*nh);

  for (int i = 0; i < nd; i++){
    double tmp = 0;
    for (int j = 0; j < nh; j++){
      if(arma::as_scalar((X.row(dinx[i])-X.row(hinx[j]))*beta)>0) tmp += 1;
    }
    out1 += (tmp/nh - eauc)*(tmp/nh - eauc);
  }
  out1 = out1/(nd*nd);

  for (int i = 0; i < nh; i++){
    double tmp = 0;
    for (int j = 0; j < nd; j++){
      if(arma::as_scalar((X.row(dinx[j])-X.row(hinx[i]))*beta)>0) tmp += 1;
    }
    out2 += (tmp/nd - eauc)*(tmp/nd - eauc);
  }
  out2 = out2/(nh*nh);

  out2 += out1;

  return out2;
}

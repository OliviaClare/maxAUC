#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>

// [[Rcpp::export]]
double eauc_l1(arma::vec beta, arma::mat X, arma::vec Y, bool silence){
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
  
  for (int i = 0; i < nd; i++){
    for (int j = 0; j < nh; j++){
      if(arma::as_scalar((X.row(dinx[i])-X.row(hinx[j]))*beta)>0) out += 1;
    }
  }
  
  out = out/(nd*nh);
  
  return out;
}

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
  
  return out;
}


// [[Rcpp::export]]
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




// [[Rcpp::export]]
double varauc_l1(arma::vec beta, arma::mat X, arma::vec Y){
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




// 
// double varauc_ustat(arma::vec beta, arma::mat X, arma::vec Y){
//   double hatauc = 0; double hatF = 0;
//   double out = 0;
//   int d = X.n_cols - 1;
//   
//   int beta_d = beta.n_elem;
//   if(beta_d!=(d+1)){
//     printf("The dimension of coefficient is incorrect.\n");
//     return 0;
//   }
//   arma::uvec dinx = find(Y==1);
//   arma::uvec hinx = find(Y==0);
//   int nd = dinx.n_elem;
//   int nh = hinx.n_elem;
//   int n = nd+nh;
//   
//   hatauc = eauc_sort(beta, X, Y, TRUE);
//   
//   hatF = n*(n-1)/(nd*nh); 
//   
//   for (int i = 0; i < n; i++){
//     double tmp = 0;
//     for (int j = 0; j < n; j++){
//       arma::vec risk = (X.row(i)-X.row(j))*beta;
//       // risk.print("risk=");
//       tmp = tmp + (arma::as_scalar(risk)>0 - hatauc)*(Y(i)>Y(j)) + (arma::as_scalar(risk)<0 - hatauc)*(Y(j)>Y(i));
//       // printf("tmp=%d", tmp);
//     }
//     out += pow((tmp/n),2);
//   }
//   
//   out = out/(n*(n-1))*hatF*hatF;
//   
//   return out;
// }


// [[Rcpp::export]]
arma::mat dtauc_opt(arma::vec beta, arma::vec beta_k, arma::mat X, arma::vec Y, double sigma0, double w, double t){
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
  out = out - w * 2 * beta_k;
  // out.print("out="); beta_k.print("beta=");
  out = out + 2*beta/(t*(1-pow(norm(beta),2)));
  return out;
}



// // [[Rcpp::export]]
// arma::mat ddtauc_sort(arma::vec beta, arma::mat X, arma::vec Y, double sigma0){
//   // sort optimized version
//   int d = X.n_cols;
//   
//   arma::mat out; out=arma::zeros(d,d);
//   arma::uvec dinx = find(Y==1);
//   arma::uvec hinx = find(Y==0);
//   arma::mat delta;
//   int nd = dinx.n_elem;
//   int nh = hinx.n_elem;
//   arma::mat xd = X.rows(dinx);
//   arma::mat xh = X.rows(hinx);
//   // xd.print("xd=");
//   
//   // calculate and sort the linear combination for each group
//   arma::vec dlc = sort(xd*beta); //default is ascending: sort_direction="ascend"
//   arma::vec hlc = sort(xh*beta);
//   // dlc.print("dlc=");
//   
//   xd = xd.rows(sort_index(xd*beta));
//   xh = xh.rows(sort_index(xh*beta));
//   
//   // xd.print("xd");
//   
//   int start_index = 0;
//   int hit = 0; int reset = 0;
//   
//   for (int i = 0; i < nd; i++){
//     if(hit==0){
//       if(dlc[i]-hlc[0]<0) continue;
//       for  (int j = 0; j <nh; j++){
//         if((dlc[i]-hlc[j] < sigma0) && (dlc[i]-hlc[j] >= 0)){
//           if(hit==0){
//             start_index = j; //starting index found
//             hit=1; // 
//           }
//           arma::mat tmp_x = xd.row(i)-xh.row(j);
//           out = out + tmp_x.t()*tmp_x;
//         }
//       }	
//     }
//     else{
//       reset = 0;
//       for  (int j = start_index; j <nh; j++){
//         if((dlc[i]-hlc[j] < sigma0) && (dlc[i]-hlc[j] >= 0)){
//           if(reset==0){
//             start_index = j; //starting index found
//             reset=1; // 
//           }
//           arma::mat tmp_x = xd.row(i)-xh.row(j);
//           out = out + tmp_x.t()*tmp_x;
//         }
//         else if(dlc[i]-hlc[j] < 0) break;
//       }
//     }
//   }
//   out = out/(nd*nh*pow(sigma0,2));
//   return out;
// }


// // [[Rcpp::export]]
// arma::vec newton_raphson_tri_opt(arma::vec beta_k, double sigma0, double tol, int iteration, double gamma, arma::mat X, arma::vec Y, double w, double t){
//   // beta_k is the coef from previous CCP outer loop
//   
//   arma::vec beta1;
//   arma::vec beta0; beta0=normalise(beta_k)*0.5;
//   arma::vec beta_half;
//   int d = beta_k.n_elem;
//   arma::mat stp;
//   beta_k = normalise(beta_k)*0.5;
//   
//   arma::mat Iddqn = arma::eye(d,d);
// 
//   arma::mat beta_store = arma::zeros(iteration,d);
//   arma::vec qn = arma::zeros(iteration);
//   for(int i=0; i < iteration; i++){
//     // printf("This is %d iteration. \n", i);
//     arma::mat ddqn = ddtauc_sort(beta0, X, Y, sigma0) + (4*beta0/pow((1-pow(norm(beta0),2)),2)*beta0.t() + Iddqn*2/(1-pow(norm(beta0),2)))/t;
//     arma::vec dqn = dtauc_opt(beta0, beta_k, X, Y, sigma0, w, t) + 2*beta0/(t*(1-pow(norm(beta0),2)));
//     
//     // stp = pinv(ddqn)*dqn; //use moore-penrose inverse for ddqn
//     
//     // arma::mat A = arma::join_rows(arma::join_vert(ddqn, arma::conv_to<arma::mat>::from(2*beta0.t())), 
//     //                               arma::join_vert(arma::conv_to<arma::mat>::from(2*beta0), arma::zeros(1)));
//     // arma::mat B = arma::join_vert(-dqn, arma::zeros(1));
//     // // ddqn.print("ddqn=");
//     // dqn.print("dqn=");
//     // beta0.print("beta0=");
//     // arma::mat weired =  (4*beta0/pow((1-pow(norm(beta0),2)),2)*beta0.t() + Iddqn*2/(1-pow(norm(beta0),2)))/t;
//     // weired.print("barrier=");
//     // A.print("A=");
//     // B.print("B=");
//     // if(ddqn.is_zero()){
//       
//       // }
//     stp = pinv(ddqn)*dqn; //use moore-penrose inverse for ddqn
// 
//     // arma::mat stp2 = pinv(A)*B;
//     // stp.print("stp=");
//     // stp2.print("stp2=");
//     // stp = stp.rows(0,d-1);
//     // stp.print("stp=");
//     stp = normalise(stp)/10;
//     
//     beta1 = beta0 - stp;
//     beta_half = beta0 - gamma*stp;
// 
//     // if(norm(beta1)>=1) {beta1.print("bigger"); beta1 = normalise(beta1)*0.8;}
//     // else beta1.print("smaller");
//     if(norm(beta1)>=1) {beta1 = normalise(beta1)*0.8;}
//     if(norm(beta_half)>=1) beta1 = normalise(beta_half)*0.8;
//     
//     // beta1.print("beta1");
//     // beta_half.print("beta_half");
//     
//     // double qn0 = eauc_sort(beta1, X, Y, TRUE);
//     // double qn1 = eauc_sort(beta_half, X, Y, TRUE);
//     double qn0 = -tauc_sort(beta1, X, Y, sigma0, TRUE) - w*norm(beta1) - log(1-pow(norm(beta0),2))/t;
//     double qn1 = -tauc_sort(beta_half, X, Y, sigma0, TRUE) - w*norm(beta_half)- log(1-pow(norm(beta0),2))/t;
//     
//     
//     if(qn0 < qn1) {
//       beta1 = beta_half; // half the step size if see no improvement
//     }
//     
//     beta_store.row(i) = beta1.t();
//     qn[i] = eauc_sort(beta1, X, Y, TRUE);	
//     // dqn = dtauc_opt(beta1, beta_k, X, Y, sigma0, w);
//     
//     // qn.print("qn=");
//     if(i>5 && accu(abs(stp))< tol*100 ){
//       //dqn.print("dqn=");
//       //th1.print("Local Maximum: ");
//       if(abs(qn[i-1]-qn[i])<tol){
//         // return i;
//         return beta1;
//       } 
//       else {
//         arma::uvec tmp = find(qn==max(qn));
//         // tmp.print("max qn index:");
//         // qn.print("qn:");
//         // th_store.print("history of theta:");
//         beta1 = beta_store.row(tmp(0)).t();
//         
//       }
//     }
//     // if(accu(abs(th1))>50) th1=arma::randu(d,1);
//     beta0 = beta1;
//   }
//   printf("Warning: Newton Raphson did not converge within %d iterations. \n", iteration);
//   arma::uvec tmp = find(qn==max(qn));
//   //tmp.print("max index=");
//   beta1 = beta_store.row(tmp(0)).t();
//   //th1.print("th1=");
//   // return iteration;
//   return beta1;
// }

arma::vec ddnormcpp(arma::vec x){
  arma::vec out;
  out = -(x*arma::normpdf(x));
  return out;
}



//update Dn
// [[Rcpp::export]]
arma::mat vn(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor){
  // beta should have anchor set to 1
  beta = beta/beta[anchor-1];
  int d = X.n_cols - 1;
  arma::mat out; out=arma::zeros(d,d); 
  arma::vec tmp;
  // arma::vec beta = add_one(theta,1,anchor);
  arma::uvec v = find(beta!=1); //index vector for non-anchor biomarkers
  arma::uvec dinx = find(Y==1);
  arma::uvec hinx = find(Y==0);
  int nd = dinx.n_elem;
  int nh = hinx.n_elem;
  int n = nd + nh;
  
  for(int i=0; i<nd; i++){
    tmp.zeros(d);
    for(int j=0; j<nh; j++){
      arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]); //row vector 1*(d+1)
      arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
      if(tmp_x.elem(v).is_zero()){
        // out += 0.5;
        continue;
      }
      tmp = tmp + tmp_x.elem(v)*(normpdf(sqrt(n)*tmp_x*beta/tmp_sig)*sqrt(n)/tmp_sig);
    }
    out = out + tmp*tmp.t();
  }
  
  for(int i=0; i<nh; i++){
    tmp.zeros();
    for(int j=0; j<nd; j++){
      arma::mat tmp_x = X.row(hinx[i])-X.row(dinx[j]); //row vector 1*(d+1)
      arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
      if(tmp_x.elem(v).is_zero()){
        // out += 0.5;
        continue;
      }
      tmp = tmp + tmp_x.elem(v)*(normpdf(sqrt(n)*tmp_x*beta/tmp_sig)*sqrt(n)/tmp_sig);
    }
    out = out + tmp*tmp.t();
  }
  
  out = out/pow(n,3);
  // out = out*n/(pow(nd,2)*pow(nh,2));
  return out;
}



// [[Rcpp::export]]
arma::mat an(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor){
  // calculates hessian
  // beta should have anchor set to 1
  beta = beta/beta[anchor-1];
  int d = X.n_cols - 1;
  arma::mat out; out=arma::zeros(d,d);
  arma::vec tmp;
  // arma::vec beta = add_one(theta,1,anchor);
  arma::uvec v = find(beta!=1); //index vector for non-anchor biomarkers
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
      arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
      // tmp_sig.print("tmp_sig");
      if(tmp_sig.is_zero()){
        continue;
      }
      mtmp = tmp_x.elem(v)*(n*ddnormcpp(sqrt(n)*tmp_x*beta/tmp_sig)/square(tmp_sig))*tmp_x.elem(v).t();
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
  eigval.elem(find(eigval>0)).zeros();
  eigval.elem(find(abs(eigval)<1e-4)).zeros();
  // eigval.print("eigval");
  // eigvec.print("eigvec");
  // // eigval.print("eigval");
  out = eigvec*diagmat(eigval)*eigvec.t();
  // out.print("out2");
  // if((-out).is_sympd()){
    // 	// printf("neg def \n");
    // }else{
      // 	// out.print("out=");
      
      // }
  
  
  return out;
}


// [[Rcpp::export]]
arma::mat dn(arma::vec beta, arma::mat X, arma::vec Y, arma::mat var, int anchor){
  // int n = X.n_rows;
  arma::mat an_inv;
  arma::mat out;
  beta = beta/beta[anchor-1];
  // beta.print("beta=");
  an_inv = pinv(an(beta,X,Y,var,anchor));
  out = an_inv*vn(beta,X,Y,var,anchor)*an_inv;
  return out;
}



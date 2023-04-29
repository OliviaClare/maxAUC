#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]


// arma::vec find_equal(arma::vec A, double b){
// 	arma::uvec idx = find(A == b); // Substitute == with >, <, >=, <=, !=
// 	arma::vec out = A.elem(idx);   // Retrieve elements from positional index
// 	return out;
// }


arma::vec add_one(arma::vec theta, int indicator){
	//append one at the end of theta
	int d = theta.n_rows;
	arma::vec out;
	if (indicator==1) out.ones(d+1);
	else out=arma::zeros(d+1);
	for (int i=0; i<d; i++){
		out[i] = theta[i];
	}
	return out;
}


// [[Rcpp::export]]
double eauc(arma::vec theta, arma::mat X, arma::vec Y){
	double out = 0;
	int d = X.n_cols - 1;
	// arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	// arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
	int theta_d = theta.n_elem;
	if(theta_d!=d){
		printf("The dimension of coefficient is incorrect.\n");
		return 0;
	}
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	//hinx.print("hinx=");
	arma::vec beta = add_one(theta,1);

	for (int i = 0; i < nd; i++){
		for (int j = 0; j < nh; j++){
			if(arma::as_scalar((X.row(dinx[i])-X.row(hinx[j]))*beta)>0) out += 1;
		}
	}

	out = out/(nd*nh);
	if(out < 0.5) out=1-out;

	return out;
}


double qncpp(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	double out = 0;
	int d = X.n_cols - 1;
	arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
	arma::vec beta = add_one(theta,1);
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	int n = nd + nh;
	//hinx.print("hinx=");

	for (int i = 0; i < nd; i++){
		for (int j = 0; j < nh; j++){
			arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]);
			//tmp_x.print("tmp_x=");
			arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
			arma::mat mtmp = sqrt(n)*tmp_x*beta;
			if(mtmp.is_zero() && tmp_sig.is_zero()) continue;
			// tmp_sig.print("tmp_sig");
			// mtmp.print("mtmp=");
			out = out + arma::normcdf(mtmp[0]/tmp_sig[0]);
		}
	}

	out = out/(nd*nh);
	return out;
}



// [[Rcpp::export]]
arma::vec dqncpp(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	int d = X.n_cols - 1;
	arma::vec out; out=arma::zeros(d);
	arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
	arma::vec beta = add_one(theta,1);
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	arma::vec mtmp;
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	int n = nd + nh;

	for (int i = 0; i < nd; i++){
		for (int j = 0; j < nh; j++){
			arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]);
			arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
			if(tmp_sig.is_zero()){
				continue;
			}
			mtmp = sqrt(n)*tmp_x*beta; //inside the arma::normpdf
			out = out + tmp_x.elem(v)*arma::normpdf(mtmp[0]/tmp_sig[0])*sqrt(n)/tmp_sig[0];
			// if(i==0){
			// 	tmp_x.print("tmp_x=");
			// 	tmp_sig.print("tmp_sig");
			// 	mtmp.print("mtmp=");
			// 	out.print("out=");
			// }
		}
		// out.print("out=");
	}
	out = out/(nd*nh);
	return out;
}



arma::vec ddnorm(arma::vec x){
	arma::vec out;
	out = -(x*arma::normpdf(x));
	return out;
}


arma::mat ddqncpp(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	int d = X.n_cols - 1;
	arma::mat out; out=arma::zeros(d,d);
	arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	arma::vec beta = add_one(theta,1);
	arma::uvec v = find(beta!=1);
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	int n = nd + nh;
	arma::mat mtmp; mtmp=arma::zeros(d,d);

	for (int i = 0; i < nd; i++){
		for (int j = 0; j < nh; j++){
			arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]);
			arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
			if(tmp_sig.is_zero()){
				continue;
			}
			//tmp_sig.print("tmp_sig:");
			mtmp = tmp_x.elem(v)*(n*ddnorm(sqrt(n)*tmp_x*beta/tmp_sig)/square(tmp_sig))*tmp_x.elem(v).t();
			out = out + mtmp;
			// out.print("out=");
		}
	}
	out = out/(nd*nh);
	// out = out/(n*(n-1));

	// //correct Hessian to semi-negative definite
	// arma::vec eigval;
	// arma::mat eigvec;
	// eig_sym(eigval, eigvec, out);
	// // eigval.print("eigval");
	// // eigvec.print("eigvec");
	// // out.print("out");
	// eigval.elem(find(eigval>0)).zeros();
	// eigval.elem(find(abs(eigval)<1e-4)).zeros();
	// // eigval.print("eigval");
	// out = eigvec*diagmat(eigval)*eigvec.t();

	// // if((-out).is_sympd()){
	// // 	// printf("neg def \n");
	// // }else{
	// // 	// out.print("out=");

	// // }

	return out;
}


// [[Rcpp::export]]
arma::mat an_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	int d = X.n_cols - 1;
	arma::mat out; out=arma::zeros(d,d);
	arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	arma::vec beta = add_one(theta,1);
	arma::uvec v = find(beta!=1);
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	int n = nd + nh;
	arma::mat mtmp; mtmp=arma::zeros(d,d);

	for (int i = 0; i < nd; i++){
		for (int j = 0; j < nh; j++){
			arma::mat tmp_x = X.row(dinx[i])-X.row(hinx[j]);
			arma::mat tmp_sig = sqrt(tmp_x.elem(v).t()*var*tmp_x.elem(v));
			if(tmp_sig.is_zero()){
				continue;
			}
			//tmp_sig.print("tmp_sig:");
			mtmp = tmp_x.elem(v)*(n*ddnorm(sqrt(n)*tmp_x*beta/tmp_sig)/square(tmp_sig))*tmp_x.elem(v).t();
			out = out + mtmp;
			// out.print("out=");
		}
	}
	// out = out/(nd*nh);
	out = out/(n*(n-1));

	// //correct Hessian to semi-negative definite
	// arma::vec eigval;
	// arma::mat eigvec;
	// eig_sym(eigval, eigvec, out);
	// // eigval.print("eigval");
	// // eigvec.print("eigvec");
	// // out.print("out");
	// eigval.elem(find(eigval>0)).zeros();
	// eigval.elem(find(abs(eigval)<1e-4)).zeros();
	// // eigval.print("eigval");
	// out = eigvec*diagmat(eigval)*eigvec.t();

	// // if((-out).is_sympd()){
	// // 	// printf("neg def \n");
	// // }else{
	// // 	// out.print("out=");

	// // }

	return out;
}



//update Dn
// [[Rcpp::export]]
arma::mat vn_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	int d = X.n_cols - 1;
	arma::mat out; out=arma::zeros(d,d);
	arma::vec tmp = arma::regspace(0,1,d-1); //index vector 1:d
	arma::uvec v = arma::conv_to<arma::uvec>::from(tmp);
	arma::vec beta = add_one(theta,1);
	arma::uvec dinx = find(Y==1);
	arma::uvec hinx = find(Y==0);
	int nd = dinx.n_elem;
	int nh = hinx.n_elem;
	int n = nd + nh;

	for(int i=0; i<nd; i++){
		tmp.zeros();
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

	// out = out/(pow(nd,2)*pow(nh,2));
	out = out/pow(n,3);
	return out;
}



// [[Rcpp::export]]
arma::mat dn_anchor(arma::vec theta, arma::mat X, arma::vec Y, arma::mat var){
	// int n = X.n_rows;
	arma::mat an_inv;
	arma::mat out;
	an_inv = pinv(an_anchor(theta,X,Y,var));
	out = an_inv*vn_anchor(theta,X,Y,var)*an_inv;
	return out;
}


// [[Rcpp::export]]
arma::mat newton_raphson_anchor(arma::vec theta_initial, double tol, int iteration, double gamma, arma::mat X, arma::vec Y, arma::mat var){
	arma::vec th0 = theta_initial; //initial value of theta, dimension d*1
	arma::vec th1;
	int d = theta_initial.n_rows;
	arma::mat stp;
	// arma::mat oned = arma::ones(d);
	// arma::vec tmp = regspace(0,1,d-1); //index vector 1:d
	// arma::uvec v = conv_to<arma::uvec>::from(tmp);
	// arma::vec theta = beta.elem(v);

	arma::mat th_store = arma::zeros(iteration,d);
	arma::vec qn = arma::zeros(iteration);
	for(int i=0; i < iteration; i++){
		arma::mat ddqn = ddqncpp(th0, X, Y, var);
		arma::vec dqn = dqncpp(th0, X, Y, var);

		if(det(ddqn)==0) stp = arma::randu(d,1);
		else	stp = pinv(ddqn)*dqn; //use moore-penrose inverse for ddqn

		th1 = th0 - stp;
		double qn0 = eauc(th1, X, Y);
		double qn1 = eauc(th0 - gamma*stp, X, Y);

		if(qn0 < qn1) {
			th1 = th0 - gamma*stp; // half the step size if see no improvement
		}

		th_store.row(i) = th1.t();
		qn[i] = eauc(th1, X, Y);
		dqn = dqncpp(th1, X, Y, var);
		//th1.print("theta.1="); stp.print("stp="); printf("qn = %f \n",qn[i]);

		if(accu(abs(dqn))<tol){
			//dqn.print("dqn=");
			//th1.print("Local Maximum: ");
			if(abs(max(qn)-qn[i])<1e-05) return(th1);
			else {
				arma::uvec tmp = find(qn==max(qn));
				// tmp.print("max qn index:");
				// qn.print("qn:");
				// th_store.print("history of theta:");
				th1 = th_store.row(tmp(0)).t();

			}
		}
		// if(accu(abs(th1))>50) th1=arma::randu(d,1);
		th0 = th1;
	}
	printf("Error: Newton Raphson did not converge within %d iterations. \n", iteration);
	arma::uvec tmp = find(qn==max(qn));
	//tmp.print("max index=");
	th1 = th_store.row(tmp[0]).t();
	//th1.print("th1=");
	return th1;
}


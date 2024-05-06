\name{aucmax}
\alias{aucmax}
\title{max AUC estimator}
\description{empirical maximization of AUC}
\usage{aucmax(X, Y,sigma1=1, sig.eps=1e-2, outer=60, eps=1e-8,
                            beta.0=NA, alpha=0.9, w=2,
                            var.method="selfinduced", tol=2, inf_iter = 10,
                            silence=FALSE)}
\arguments{
  \item{X}{matrix of biomarker values, with each row
  containing multiple markers from an individual.}
  \item{Y}{vector of disease status}
  \item{sigma1}{initial bandwidth for smoothing empirical AUC}
  \item{sig.eps}{The absolute convergence tolerance for sigma (bandwidth)}
  \item{outer}{maximum iteration for outer loop in point estimation}
  \item{eps}{The absolute convergence tolerance for AUC}
  \item{beta.0}{initial coefficient estimate}
  \item{alpha}{reduction factor in the sequence of approximation parameter sigma.}
  \item{w}{weight for l2 norm of combination coefficient in the objective
  function (w>sqrt(2) guarantees sound asymptotic properties).}
  \item{var.method}{Whether to compute variance estimate, can be NA or "selfinduced". Default value is "selfinduced".}
  \item{tol}{The absolute convergence tolerance for covariance matrix}
  \item{inf_iter}{maximum iteration in variance estimation}
  \item{silence}{whether to print intermediate results}

}
\value{
  \item{eAUC}{empirical AUC.}
  \item{beta}{estimated combination coefficient, with unity l2 norm.}
  \item{time.elapsed}{time spent for maximization}
  \item{var.auc}{variance estimate for empirical AUC}
  \item{bias_correction}{bias corrected and logit transformed bias corrected predictive and ideal AUC and its confidence interval}
  \item{cov}{covariance estimate for coefficent}
}
\examples{
## simulate 3 biomarkers
dat.train <- sim4(100, markers=3, scenario = 1, pd = 0.3, design="cohort")
Y = dat.train$y
X = as.matrix(dat.train[,-1])
aucmax(X, Y)
}
\author{Yuxuan Chen}
\note{Non-commercial use only. No warranty is implied. Please report bugs
  to Yuxuan Chen (yuxuan.chen@emory.edu).}

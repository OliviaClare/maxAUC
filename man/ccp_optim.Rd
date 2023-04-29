\name{auccp}
\alias{auccp}
\title{max AUC estimator}
\description{empirical maximization of AUC}
\usage{auccp(X, Y, sigma1=0.5, iter=20, outer=60, eps=1e-5, sig.eps=1e-2,
                         beta.0=NA, tol=1e-2,
                         alpha=0.8, alphat=0.8, w=2, t=1/sigma1, silence=FALSE)}
\arguments{
  \item{X}{matrix of biomarker values, with each row
  containing multiple markers from an individual.}
  \item{Y}{vector of disease status}
  \item{sigma1}{initial bandwidth for smoothed AUC}
  \item{iter}{maximum iteration for vaiance estimation}
  \item{outer}{maximum iteration for outer loop in point estimation}
  \item{eps}{The absolute convergence tolerance for AUC}
  \item{sig.eps}{The absolute convergence tolerance for sigma}
  \item{beta.0}{initial coefficient estimate.}
  \item{tol}{The absolute convergence tolerance for covariance matrix}
  \item{alpha}{reduction factor in the sequence of approximation parameter sigma.}
  \item{alphat}{reduction factor in the sequence of approximation parameter t.}
  \item{w}{weight for l2 norm of combination coefficient in the objective
  function (w>sqrt(2) guarantees sound asymptotic properties).}
  \item{t}{initial value for t}
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
auccp(X, Y)
}
\author{Yuxuan Chen}
\note{Non-commercial use only. No warranty is implied. Please report bugs
  to Yuxuan Chen (yuxuan.chen@emory.edu).}

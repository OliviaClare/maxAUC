\name{smrce}
\alias{smrce}
\title{max AUC estimator}
\description{empirical maximization of AUC implemented in Zhang et al. (2018)}
\usage{smrce(X,Y, beta_initial, tol,iteration, silence = TRUE)}
\arguments{
  \item{X}{matrix of biomarker values, with each row
  containing multiple markers from an individual.}
  \item{Y}{vector of disease status}
  \item{beta_initial}{initial coefficient estimate.}
  \item{tol}{The absolute convergence tolerance}
  \item{iteration}{maximum iteration for estimation}
  \item{silence}{whether to print intermediate results}

}
\value{
  \item{empirical_AUC}{maximized empirical AUC.}
  \item{beta}{estimated combination coefficient, with unity l2 norm.}
  \item{time.elapsed}{time spent for maximization.}
  \item{var.auc}{variance estimate for empirical AUC.}
  \item{cov}{covariance estimate for coefficent.}
  \item{iter}{number of iterations ran.}
  \item{converged}{status of convergence, 1 for converged.}
}
\references{Zhang, J., Jin, Z., Shao, Y., & Ying, Z. (2018). Statistical inference on transformation models: a self-induced smoothing approach. Journal of Nonparametric Statistics, 30(2), 308-331.}
\examples{
## simulate 3 biomarkers
dat.train <- sim4(100, markers=3, scenario = 1, pd = 0.3, design="cohort")
Y = dat.train$y
X = as.matrix(dat.train[,-1])
b <- glm(Y~X, family = binomial())
b.beta <- b$coefficients[-1]/norm(b$coefficients[-1], type="2")

smrce(X, Y, b.beta,1e-5,20,silence=FALSE)
}
\author{Yuxuan Chen}
\note{Non-commercial use only. No warranty is implied. Please report bugs
  to Yuxuan Chen (yuxuan.chen@emory.edu).}

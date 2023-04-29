\name{maxAUC-package}
\alias{maxAUC-package}
\alias{maxAUC}
\docType{package}
\title{
\packageTitle{maxAUC}
}
\description{
\packageDescription{maxAUC}
}
\details{
Provides fast maximization of empirical AUC for biomarker combination. Performance estimation
including bias correction and interval estimation are provided as well.
}
\author{
\packageAuthor{maxAUC}

Maintainer: \packageMaintainer{maxAUC}
}
\references{
~~ Literature or other references for background information ~~
}

\keyword{AUC, biomarker combination, bias correction, interval estimation}
\examples{
dat.train <- sim4(100, markers=3, scenario = 1, pd = 0.3, design="cohort")
Y = dat.train$y
X = as.matrix(dat.train[,-1])
auccp(X, Y)
}

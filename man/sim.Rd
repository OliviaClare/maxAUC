\name{sim4}
\alias{sim4}
\title{data simulation}
\description{empirical maximization of AUC}
\usage{sim4(sample.size=100, scenario=1, markers=3, pd=0.3, design="cohort")}
\arguments{
  \item{sample.size}{sample size}
  \item{scenario}{scenario range from 1 to 4.}
  \item{markers}{number of markers, can be 3 or 6}
  \item{pd}{prevalence of disease}
  \item{design}{"cohort" for cohort design, "casecontrol" for case control design}

}
\value{
  \item{res}{a data frame of biomarker values X and disease status Y}
}
\examples{
## simulate 3 biomarkers for 100 cases and 100 controls
dat.train <- sim4(100, markers=3, scenario = 1, pd = 0.3, design="cohort")
}
\author{Yuxuan Chen}
\note{Non-commercial use only. No warranty is implied. Please report bugs
  to Yuxuan Chen (yuxuan.chen@emory.edu).}

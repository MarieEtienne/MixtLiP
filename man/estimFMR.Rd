% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/FullMixtureRegressionFunctions.R
\name{estimFMR}
\alias{estimFMR}
\title{estimFMR runs the EM algorithm to estimate the FMR model paraters}
\usage{
estimFMR(data, K, initParameters = NULL, nIter = 10, garde = TRUE,
  file = "parameters")
}
\arguments{
\item{data}{the data (d columns for X and the last on e for Y)}

\item{K}{the number of cluster}

\item{initParameters}{optionnal set of parameters to  start EM algorithm}
}
\value{
the parameters after running the EM algorithm
}
\description{
estimFMR runs the EM algorithm to estimate the FMR model paraters
}


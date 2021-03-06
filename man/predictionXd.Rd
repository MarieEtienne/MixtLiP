% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/EstimationPredictionFunctions.R
\name{predictionXd}
\alias{predictionXd}
\title{predictionXd predicts Xd given the estimated parameters of the mixture and X1:(d-1)}
\usage{
predictionXd(mixtModel, covariate, Y = NULL)
}
\arguments{
\item{mixtModel}{an object of class MixmodResults containing the mixmodel used for the prediction}

\item{covariate}{the vector (or matrix (n,d-1) ) X1:(d-1) of the covariates to predict Xd}
}
\value{
a numeric vector Xd size(n)
}
\description{
predictionXd predicts Xd given the estimated parameters of the mixture and X1:(d-1)
}


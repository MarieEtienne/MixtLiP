% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BootstrapFunctions.R
\name{booStrapStabilityforK}
\alias{booStrapStabilityforK}
\title{booStrapStabilityforK chechs stability}
\usage{
booStrapStabilityforK(data, Kmax = 5, nBoot = 10, gBoot = FALSE,
  group = NULL, crit = "ICL", fileName = "KStat", model = "all",
  resDir = "", mc.cores = 2)
}
\arguments{
\item{data}{the data to  boostrap from}

\item{Kmax}{the maximum number of components}

\item{nBoot}{the number of boostrap sampling}

\item{gBoot}{logical if TRUE a parametric boostrap who preserves group structure, and in this case argument group as to be specified, 
if FALSE a force brut boostrap (choose n among n)}

\item{group}{a vector specifying a group for each individual, requires for parametric bootstrap}

\item{crit}{the criterion to be maximixed for the selection of the number of components}

\item{fileName}{the name uses to store the results}

\item{model}{the kind of model to be tested}

\item{resDir}{the optionnal directory name to store the results}

\item{mc.cores}{the number of cores to be used default is 2}
}
\value{
a nBoot length vector containing the number of groups selected
}
\description{
booStrapStabilityforK chechs stability
}


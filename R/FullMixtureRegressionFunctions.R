# initEMFullMixtureRegression
#' initEMFullMixtureRegression uses Rmixmod  and regression to provide good initialisations for the FMR EM algorithm.
#' @param data the data to be fitted
#' @param K the number of cluster
#' @import Rmixmod
#' @export
#' @return a list of parameters with
#'  pik the proportion of the mixture,
#'  sigmaX a list of K covariance matrix,
#'  muX  a list of K d dimensionnal expectation and
#'  theta a list of K regression parameters
#'  sigmaY a list of K residual variance for the regression
initEMFullMixtureRegression <- function(data, K){
  data <- as.data.frame(data)
  d <- ncol(data)-1
  names(data) <- c(paste0('X', 1:d),'Y')

  X <- data[,1:d]
  Y <- data[,d+1]
  n <- nrow(data)
  X.clustering <- mixmodCluster(data = as.data.frame(X), nbCluster = K,
                                model =  mixmodGaussianModel(family = "general", listModels = 'Gaussian_pk_Lk_Ck',
                                                             free.proportions = TRUE, equal.proportions = FALSE))

  sigmaX <- X.clustering@bestResult@parameters@variance
  muX <- lapply(1:K,
                function(k){X.clustering@bestResult@parameters@mean[k,]})
  pik <- X.clustering@bestResult@parameters@proportions
  tau <- X.clustering@bestResult@proba

  form <- paste0('Y~', names(data)[1], '+', names(data)[2])
  beta <- lapply(1:K, function(k){
    dataprov <- data.frame(data, w=tau[,k])
    lmProv <- lm(form,  weights =w, data=dataprov)
    return(coef(lmProv))
  }
  )

  sigmaY<- lapply(1:K, function(k){
    dataprov <- data.frame(data, w=tau[,k])
    lmProv <- lm(form,  weights =w, data=dataprov)
    return(summary(lmProv)$sigma^2)
  })

  return(list(muX=muX, sigmaX=sigmaX, beta=beta, sigmaY = sigmaY, pik=pik))
}

# drawFMRParam
#' drawFMRParam draws a set of parameter for the FMR model
#' @param K the number of cluster
#' @param d the dimension for X
#' @param parameters the parameter that is a list with component
#' pik the mixture proportion
#' foreach k in 1,\ldots, K mu_k, \Sigma_k, \beta_k and \sigma_k^Y
#' if parameters == NULL, then the set of parameters is drawn
#' @export
#' @return a list with Y, X the observations, Z and parameters the list of parameters
drawFMRParam <- function(K, d, p=3){
  mu  <-
    lapply(1:K, function(i) {
      round(matrix(3 * rnorm(d), ncol = 1),1)
    }) ## mean values
  ## precision matrix, randomly drawn from Wishart distribution
  Gamma <- lapply(1:K, function(i) {
    mat <- round(matrix(rnorm(p * d, mean = 0, sd = 1), nrow = p, ncol = d),1) ## wishart Wp( Id, d),
    t(mat) %*% mat
  })
  GammaInv <- lapply(Gamma, function(d) {
    solve(d)
  })
  pik <- runif(K)
  pik <- round(pik / sum(pik), 2)

  betak <- lapply(1:K,function(i){2*rnorm(d+1)})
  sigma2k <- lapply(1:K,function(i){rnorm(1)^2})

  parameters <- list(muX=mu, sigmaX=Gamma, beta=betak, sigmaY = sigma2k, pik=pik, sigmaXInv=GammaInv )
  return(parameters)
}

# simFMR
#' ssimFMR draws sample from model FMR in the paper. That is a mixiture where given component k
#' X_i\vert Zik=1 is a normal distribution of mean mu_k and variance \Sigma_k and given X_i and Z_ik Y is normal
#' distribution with mean t(X_i) \beta_k and variance \sigma^Y_k
#' @param n the number of observations
#' @param K the number of cluster
#' @param d the dimension for X
#' @param parameters the parameter for the FMR model, that is a list with component
#' pik the mixture proportion
#' foreach k in 1,\ldots, K muX_k, \SigmaX_k, \beta_k and \sigmaY_k
#' if parameters == NULL, then the set of parameters is drawn
#' @import mvtnorm
#' @export
#' @return a list with Y, X the observations, Z and parameters the list of parameters
simFMR <- function(n,K, parameters=NULL,seed=NULL, d=2){
  p <- 10 # for the whishart distribution
  if(!is.null(seed)){
    set.seed(seed)
  }
  if(is.null(parameters)){
    parameters <- drawFMRParam(K = K, d = d)
  }

  if(! ('GammaInv'%in% names(parameters))){
    parameters$sigmaXInv <- lapply(1:K, function(k){solve(parameters$sigmaX[[k]])})
  }

  Z <- sample(1:K, size = n, replace = T, prob = parameters$pik) ## clussters affectation
  X <- t(sapply(1:n, function(i){rmvnorm(1, mean=parameters$muX[[Z[i]]], sigma = parameters$sigmaX[[Z[i]]])}))
  Y <- sapply(1:n, function(i){rnorm(1, mean=parameters$beta[[Z[i]]][1]+
                                       sum(parameters$beta[[Z[i]]][2:(d+1)]*X[i,]),
                                     sd = sqrt(parameters$sigmaY[[Z[i] ]])
  )})

  return(list(Y=Y, X=X, Z=Z, parameters=parameters))
}

# updateFMRParameters
#' updateParameters runs one iteration of the EM algorithm
#' @param K the number of cluster
#' @param data the data (d columns for X and the last on e for Y)
#' @param currentParameters the current state of parameters (theta')
#' @export
#' @return the result of the optimisation of Q(theta, theta')
updateFMRParameters <- function(data, K, currentParameters){
  names(data) <- c(paste0('X', 1:d),'Y')
  X <- as.data.frame(data[,1:d])
  Y <- data[,d+1]

  tau <- sapply(1:K, function(k){
    dnorm(Y, mean= currentParameters$beta[[k]][1]+ as.matrix(X, ncol=d)%*%matrix(currentParameters$beta[[k]][2:(d+1)],ncol=1), sd=sqrt(currentParameters$sigmaY[[k]]))*
      dmvnorm(data[,1:d],mean=currentParameters$muX[[k]], sigma=currentParameters$sigmaX[[k]])*currentParameters$pik[[k]]})
  tau = tau / rowSums(tau)
  tauPlus <- colSums(tau)

  muX <- lapply(1:K,function(k){colSums(data[,1:d]*tau[,k])/tauPlus[k]})

  sigmaX <- lapply(1:K,function(k){
    Wk <- sweep(as.matrix(X), MARGIN = 2, STATS = muX[[k]], FUN = '-')
    Reduce(f = '+', x = lapply(1:n, function(i){
      matrix(as.numeric(Wk[i,]), ncol=1) %*% matrix(as.numeric(Wk[i,]), nrow=1)*tau[i,k]
    }))/tauPlus[k]
  })


  form <- paste0('Y~', names(data)[1], '+', names(data)[2])
  beta <- lapply(1:K, function(k){
    dataprov <- data.frame(data, w=round(tau[,k],5))
    lmProv <- lm(form,  weights =w, data=dataprov)
    return(coef(lmProv))
  }
  )

  sigmaY<- lapply(1:K, function(k){
    dataprov <- data.frame(data, w=round(tau[,k],5))
    lmProv <- lm(form,  weights =w, data=dataprov)
    return(summary(lmProv)$sigma^2)
  })

  pik <- tauPlus/sum(tauPlus)
  return(list(muX=muX, sigmaX=sigmaX, beta=beta, sigmaY=sigmaY, pik=pik, tau=tau))
}


# estimFMR
#' estimFMR runs the EM algorithm to estimate the FMR model paraters
#' @param K the number of cluster
#' @param data the data (d columns for X and the last on e for Y)
#' @param initParameters optionnal set of parameters to  start EM algorithm
#' @export
#' @return the parameters after running the EM algorithm
estimFMR <- function(data, K, initParameters=NULL, nIter=10, garde=TRUE, file='parameters'){
  d <- dim(data)[2]-1
  muGarde <- matrix(NA, ncol=K*d, nrow=nIter)
  if(is.null(initParameters)){
    currentParameters <- initEMFullMixtureRegression(data = data, K=K)
  }else
  {
    currentParameters <- initParameters
  }
  if(garde){
    histParameters=list(currentParameters)
  }
  for( i in 2:nIter){
    currentParameters <- updateFMRParameters(data, K, currentParameters)
    if(garde){
      histParameters[[i]] <- currentParameters
    }
  }

  if(!garde){
    return(currentParameters)
  }else{

    return(histParameters)
  }

}

# logLikeFMR
#' logLikeFMR computes the log likelihood of the FMR model
#' @param K the number of cluster
#' @param data the data (d columns for X and the last on e for Y)
#' @param currentParameters a list containing theta
#' @export
#' @return the loglikelihood l(\theta; Y)
logLikeFMR <- function(data, K, currentParameters){
  d <- dim(data)[2]-1
  n <- data(data)[1]
  X <- as.data.frame(data[,1:d])
  Y <- data[,d]
  tau <- sapply(1:K, function(k){
    dnorm(Y, mean= currentParameters$beta[[k]][1]+ as.matrix(X, ncol=d)%*%matrix(currentParameters$beta[[k]][2:(d+1)],ncol=1), sd=sqrt(currentParameters$sigmaY[[k]]))*
      dmvnorm(X,mean=currentParameters$muX[[k]], sigma=currentParameters$sigmaX[[k]])*currentParameters$pik[[k]]
    })
  return( sum(log(rowSums(tau))) )
}

# BICFMR
#' BICFMR computes the BIC of the FMR model
#' @param data the data (d columns for X and the last on e for Y)
#' @param hatTheta a list containing  the maximum likelihood estimation
#' @export
#' @return the BIC criterio  -2 l(\theta; Y) +2 log (n) * (K-1 + (d+1) K + d K + d (d+1) /2 *K )
BICFMR <- function(data,  hatTheta){
 n <- nrow(data)
 d <- dim(data)[2]-1
 K <- length(hatTheta$muX)
 return( -2*logLikeFMR(data, K, hatTheta) +2 * log(n)*(K-1 + (d+1) * K + d * K + d* (d+1) /2 *K ) )
}

# predictFMR
#' predictFMR computes the prediction for the given covariates data with  FMR model and parameter theta
#' @param covariates a matrix or data.frame with the d columns of X
#' @param  theta a list containing  the parameter value to be used for prediction
#' @export
#' @return a vector of length nrow(covariates) containing the prediction i.e \hat{Y^{(k)}}_i \P(Z_i=k\vert X_i)
predFMR <- function(covariates,theta){
  n <- nrow(covariates)
  d <- dim(covariates)[2]
  K <- length(hatTheta$muX)
  weight <- sapply(1:K, function(k){
    dmvnorm(covariates,mean=theta$muX[[k]], sigma=theta$sigmaX[[k]])*theta$pik[[k]]
  })
  Yhat <- sapply(1:K, function(k){
    theta$beta[[k]][1]+ as.matrix(covariates, ncol=d)%*%matrix(theta$beta[[k]][2:(d+1)],ncol=1)})
  pred <- rowSums(weight*Yhat)
  return( pred)
}


# estimateMixtureModel
#' estimateMixtureModel uses Rmixmod To estimate a d dimensionnal mixture model and save the results with the appropriate fileName
#' @param data the data to be fitted
#' @param Kmax the maximum number of components
#' @param crit the criterion to be maximixed for the selection of the number of components
#' @param fileName the name uses to store the results
#' @param model the kind of model to be tested
#' @param ResDir the optionnal directory name to store the results
#' @import Rmixmod
#' @export
#' @return a Rmixmod results
estimateMixtureModel <- function(data, Kmax=5, crit='ICL', fileName='ClusterData', model='all', resDir='') {
  Y.clustering <-   mixmodCluster(data = as.data.frame(data),
                                  nbCluster = 1:Kmax,
                                  model = mixmodGaussianModel(family = model),
                                  criterion = crit)
  save('Y.clustering', file = file.path(resDir,paste0(fileName,crit,'.Rd')))
  return(Y.clustering)
}



# predictionXd
#' predictionXd predicts Xd given the estimated parameters of the mixture and X1:(d-1)
#' @param mixtModel an object of class MixmodResults containing the mixmodel used for the prediction
#' @param covariate the vector (or matrix (n,d-1) ) X1:(d-1) of the covariates to predict Xd
#' @import Rmixmod
#' @export
#' @return a numeric vector Xd size(n)

predictionXd <- function(mixtModel, covariate, Y=NULL) {
  K     <-   mixtModel@nbCluster
  Sigma <-
    mixtModel@parameters@variance ## a list with K component, each component a (d,d) matrix
  mu    <-
    mixtModel@parameters@mean ## a matrix with K rows and d columns
  pi    <- mixtModel@parameters@proportions ## a vector
  d <- length(mu[1,])

  ## prediction is of theform Y = \sum_k Xdes \theta^k pik
  ## Xdes <- (1 covariates)
  if (!is.matrix(covariate) ) {
    if(is.data.frame(covariate)){
      covariate <- as.matrix(covariate)
    }else{
    covariate <- matrix(as.numeric(covariate), ncol = d - 1)
    }
  }
  n <- nrow(covariate)
  Xdes <- matrix(c(rep(1,n),covariate), ncol = d)
  sigmaMoins1d <-
    lapply(1:K, function(k) {
      print(k)
      solve(Sigma[[k]][1:(d - 1), 1:(d - 1)])
    })
  sigmaMoins1sigmad <- lapply(1:K, function(k) {
    Sigma[[k]][d, 1:(d - 1)] %*% sigmaMoins1d[[k]]
  })
  pred <- sapply(1:K, function(k) {
    theta <-
      c(mu[k,d] - sigmaMoins1sigmad[[k]] %*% matrix(mu[1:(d - 1)], ncol = 1),  sigmaMoins1sigmad[[k]])
    Xdes %*% theta
  })
  tauk <-
    sapply(1:K, function(k) {
      dmvnorm(x = covariate, mean = mu[k,1:(d - 1)],sigma =  Sigma[[k]][1:(d -
                                                                             1), 1:(d - 1)]) * pi[k]
    })
  if (is.matrix(tauk)) {
    tauk <- sweep(tauk, MARGIN=1,STATS=rowSums(tauk), FUN = '/')
  }else{
    tauk <- tauk / sum(tauk)
  }
  if(is.matrix(tauk)){
  Y.pred <- rowSums(pred * tauk)
  }else{
    Y.pred <- sum(pred * tauk)
  }
  if(is.null(Y)){
    RMSE <- NULL
  } else{
    RMSE <- sqrt(mean((Y-Y.pred)^2))
  }
  return(list(Y.pred=Y.pred, tauk=tauk, RMSE=RMSE))
}


# estimateXonlyMixtureModel
#' estimateMixtureModel uses Rmixmod To estimate a d dimensionnal mixture model and save the results with the appropriate fileName
#' @param data the data to be fitted, the slast column being the data to be predicted
#' @param Kmax the maximum number of components
#' @param crit the criterion to be maximixed for the selection of the number of components
#' @param fileName the name uses to store the results
#' @param model the kind of model to be tested
#' @param ResDir the optionnal directory name to store the results
#' @import Rmixmod
#' @export
#' @return a list of the BestModel for clustering x and the regresssion coef
estimateXonlyMixtureModel <- function(data, Kmax=5, crit='ICL', fileName='ClusterXData', model='all',
                                      resDir='') {

  d <- ncol(data)
  Xonly.clustering <-   mixmodCluster(data = as.data.frame(data[, 1:(d-1)]),
                                  nbCluster = 1:Kmax,
                                  model = mixmodGaussianModel(family = model),
                                  criterion = crit)
  save('Xonly.clustering', file = file.path(resDir,paste0(fileName,crit,'.Rd')))
  bestModel <- Xonly.clustering['bestResult']
  K <- bestModel@nbCluster
  regCoef <- lapply(1:K, function(i){
    dataWithinClusterI <- data[which(bestModel@partition==i),]
    names(dataWithinClusterI)<- c('x1', 'x2', 'y')
    return(coef(lm(y~x1+x2, data=dataWithinClusterI)))
  })
  return(list(bestModel=bestModel, regCoef=regCoef))
}


# predictionXonly
#' predictionXonly predicts Xd given the estimated parameters of the mixture and the regression coef within component X1:(d-1)
#' @param XonlyModel a list containing bestModel, the mixture model to create the group and regCoef the regression coef within group
#' @param covariate the vector (or matrix (n,d-1) ) X1:(d-1) of the covariates to predict Xd
#' @param Y an optionnel vector with length equals to the number of rows for Xd. If given, the rmse is computed.
#' @import Rmixmod
#' @export
#' @return a numeric vector Xd size(n)

predictionXonly <- function(XonlyModel, covariate, Y=NULL) {
  K     <-   XonlyModel$bestModel@nbCluster
  Sigma <-
    XonlyModel$bestModel@parameters@variance ## a list with K component, each component a (d,d) matrix
  mu    <-
    XonlyModel$bestModel@parameters@mean ## a matrix with K rows and d columns
  pi    <- XonlyModel$bestModel@parameters@proportions ## a vector
  d <- length(mu[1,])+1

  ## prediction is of theform Y = \sum_k Xdes \theta^k pik
  ## Xdes <- (1 covariates)
  if (!is.matrix(covariate) ) {
    if(is.data.frame(covariate)){
      covariate <- as.matrix(covariate)
    }else{
      covariate <- matrix(as.numeric(covariate), ncol = d - 1)
    }
  }
  n <- nrow(covariate)
  Xdes <- matrix(c(rep(1,n),covariate), ncol = d)

  pred <- sapply(1:K, function(k) {
    Ypred <- XonlyModel$regCoef[[k]][1] + XonlyModel$regCoef[[k]][2]* covariate[,1] + XonlyModel$regCoef[[k]][3] * covariate[,2]
    return(Ypred)
})

  tauk <-
    sapply(1:K, function(k) {
      dmvnorm(x = covariate, mean = mu[k,1:(d - 1)],sigma =  Sigma[[k]]) * pi[k]
    })
  if (is.matrix(tauk)) {
    tauk <- sweep(tauk, MARGIN=1,STATS=rowSums(tauk), FUN = '/')
  }else{
    tauk <- tauk / sum(tauk)
  }
  if(is.matrix(tauk)){
    Y.pred <- rowSums(pred * tauk)
  }else{
    Y.pred <- sum(pred * tauk)
  }
  if(is.null(Y)){
    RMSE <- NULL
  } else{
    RMSE <- sqrt(mean((Y-Y.pred)^2))
  }
  return(list(Y.pred=Y.pred, tauk=tauk, RMSE=RMSE))
}


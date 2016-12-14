# booStrapStabilityforK
#' booStrapStabilityforK chechs stability
#' @param data the data to  boostrap from
#' @param Kmax the maximum number of components
#' @param nBoot the number of boostrap sampling
#' @param gBoot logical if TRUE a parametric boostrap who preserves group structure, and in this case argument group as to be specified, 
#' if FALSE a force brut boostrap (choose n among n)
#' @param group a vector specifying a group for each individual, requires for parametric bootstrap
#' @param crit the criterion to be maximixed for the selection of the number of components
#' @param fileName the name uses to store the results
#' @param model the kind of model to be tested 
#' @param resDir the optionnal directory name to store the results
#' @param mc.cores the number of cores to be used default is 2
#' @import Rmixmod, parallel
#' @export
#' @return a nBoot length vector containing the number of groups selected
booStrapStabilityforK <- function(data, Kmax=5, nBoot=10,  gBoot=FALSE, group=NULL, crit='ICL', fileName='KStat', 
                                  model='all', resDir='', mc.cores=2) {
  print(nBoot)
  if(is.null(group) & gBoot){
    cat('the group information has to be specified for grouped bootstrap.')
    return(-1)
  }
  
  n<- nrow(data)
ind <- 1:n
Klist <- mclapply(1:nBoot, function(i){
  cat(i,' ')
  print(i)
  if(gBoot){
    ind <- as.numeric(unlist(by(ind, group, function(d){ 
      sample(as.character(d), size=length(d), replace=TRUE)})))
  } else{
    ind <- sample(1:n, size=n, replace = T)  
  }
  Y <- data[ind,]
  Y.clustering <-   mixmodCluster(data = as.data.frame(data), 
                                  nbCluster = 1:Kmax, 
                                  model = mixmodGaussianModel(family = model),
                                  criterion = crit)
  K <- Y.clustering["bestResult"]@nbCluster
  mod <- Y.clustering['bestResult']@model
write(paste(K, mod, sep='\t'), file = file.path(resDir, paste0(fileName, crit, '.txt')), append = TRUE)
  return(K)
  },     mc.cores=mc.cores )
  
  return(unlist(Klist))
}
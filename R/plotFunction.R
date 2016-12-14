plot4dimData <- function(Y, Z = NULL, names=NULL, cex=1) {
  if (is.null(Z)) {
    Z = rep(1, nrow(Y))
  }
  if(is.null(names)){
    names<- paste('Var', 1:ncol(Y))
  }
  par(mfcol = c(3,3))
  plot(Y[,1], Y[,2], col = Z, xlab='',ylab=names[2], cex=cex)
  plot(Y[,1], Y[,3], col = Z, xlab='',ylab=names[3], cex=cex)
  plot(Y[,1], Y[,4], col = Z, xlab=names[1],ylab=names[4], cex=cex)
  plot(
    x = 1, axes = F, col = 0, xlab = '', ylab = '', cex=cex
  )
  plot(Y[,2], Y[,3], col = Z, xlab='', ylab='', cex=cex)
  plot(Y[,2], Y[,4], col = Z, xlab=names[2], ylab='', cex=cex)
  plot(
    x = 1, axes = F, col = 0, xlab = '', ylab = '', cex=cex
  )
  
  plot(
    x = 1, axes = F, col = 0, xlab = '', ylab = '', cex=cex
  )
  plot(Y[,3], Y[,4], col = Z, xlab=names[3], ylab='', cex=cex)
}

plot3dimData <- function(Y, Z = NULL, names=NULL, cex=1, pch=1) {
  if (is.null(Z)) {
    Z = rep(1, nrow(Y))
  }
  if(is.null(names)){
    names<- paste('Var', 1:ncol(Y))
  }
  par(mfcol = c(2,2))
  plot(Y[,1], Y[,2], col = Z, xlab='',ylab=names[2], cex=cex, pch=pch)
  plot(Y[,1], Y[,3], col = Z,ylab=names[3], xlab=names[1], cex=cex, pch=pch)
  plot(
    x = 1, axes = F, col = 0, xlab = '', ylab = '', cex=cex, pch=pch
  )
  plot(Y[,2], Y[,3], col = Z,  ylab='', xlab=names[2], cex=cex, pch=pch)
  }



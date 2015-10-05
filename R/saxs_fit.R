#' fGauss(x, p0, pf, pc, pw)
#' @description calculate Gaussian
#' @author K. Juraic
#' @param x x coord for which Gauss is calculated
#' @param p0 shift factor
#' @param pf scaling factor
#' @param pc position of gaussian peak
#' @param pw width of gaussian peak
#' @return y = Gaussina value at x
#' @examples
#' fGauss(.2,1,10,5,2)
fGauss <- function(x, p0, pf, pc, pw)
{
  out <- p0 + pf * exp(-((x-pc)^2/(2*pw^2)))
  out
}

#' fitGauss(x, y, subset)
#' @description fit Gaussian to x,y data with automatic determination of init parameters
#'              fit by using nlsLM from minpac.lm package
#' @author K. Juraic
#' @param x array of independent variable
#' @param y variable of dependent variable
#' @param subset TRUE/FALSE array which x,y pairs to include in fit
#' @return fit fited parameters
#' @examples
#' fGauss(.2,1,10,5,2)
fitGauss <- function(x, y, subset = rep(TRUE,times = length(x))){
  pInit <- list(p0=min(y) ,pf=max(y), pc=x[which.max(y)], pw=abs(max(x)-min(x))/2)
  #print(pInit)
  fit <- nlsLM( y~fGauss(x,p0,pf,pc,pw),  start = pInit, subset=subset)
  print(fit$info)
  return(fit)
}

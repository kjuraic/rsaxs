

#' Shultz-Zimm distribution
#' @author K. Juraic
#' @description Calculate Shultz Size distribution density
#'              For particle size distribution is better to use Shultz-Zihm distribution
#'              then normal size distribution because its possible to avoid negative density values
#'              SOURCE: M. khaneft, phd thesis: polymers in aligned carbon nanotube arrays (2013)
#' @param r x value for which calculation is done
#' @param Rmean distribution paramaterer (average)
#' @param Rsd distribution parameter (standard deviation)
#' @return distribution density for given x
#' @examples
#'        dSZ(2.0, 1.0, 0.3)
dSZ <- function(r, Rmean, Rsd){
  z <- (Rmean / Rsd) ^ 2 - 1
  z1 <- z + 1
  out <- (z1/Rmean) ^ z1 * r ^ z * exp(-z1 / Rmean * r) / gamma(z1)
  return(out)
}

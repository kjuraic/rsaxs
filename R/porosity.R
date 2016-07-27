#' Porosity estimation from crtical angle for total external reflectio
#' @author K. Juraic
#' @description Porosity estimation from crtical angle for total external reflectio
#' @param ac_exp experimental critical angle for total external reflection
#' @param ac_bulk theoretical bulk critical angle for ideal material
#' @return value for porosity [0-1]
#' @examples
#'  # Porosity calculation for amorphous silicon
#'  porosity(.18, .223)
porosity <- function(ac_exp, ac_bulk){
  poros <- 1-(ac_exp / ac_bulk)^2
  cat("Porosity =", poros)
  return(poros)
}



#' Cylinder Form Factor (IsGISAXS)
#' @author K. Juraic
#' @description calculate form factor for cylinder (radius R and height H) for
#'              wave vector Q(qx, qy, qz). Cylinder Z axis is normal to sample 
#'              surface.
#' @references IsGISAXS manual
#' @param qx x component of wave vector Q
#' @param qy y component of wave vector Q
#' @param qz z component of wave vector Q
#' @param R cylinder radius in nm
#' @param H cylinder height in nm
#' @return Form Factor value
#' @examples
#'    \dontrun{ff_cylinder_isgisaxs(.1, 1., .5, 10, 100)}
ff_cylinder_isgisaxs <- function(qx, qy, qz, R, H){
  qp <- sqrt(qx^2 + qy^2)
  V <- pi * R^2 * H
  bes1_par <- qp * R
  exp_par <- qz * H / 2.
  bes1 <- besselJ(bes1_par, 1) / bes1_par
  sinc <- sin(exp_par) / exp_par
  expi <- exp(1i * exp_par)
  ff <- 2 * V * bes1 * sinc * expi
  return(ff) 
}




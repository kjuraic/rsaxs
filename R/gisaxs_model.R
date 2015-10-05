#' Yoneda(tthf, af, ai, delta, beta, pI0, pF, pE, lambda = .154)
#' @description calculate intensity arround Yoneda maksimum
#' @author K. Juraic
#' @param tthf scattering angle out of specular plane
#' @param af scattering angle in specular plane
#' @param ai angle of incidence
#' @param beta index of refractong (Re part)
#' @param delta index of refraction (Im part)
#' @param surfRms surface roughness in [nm]
#' @param pI0 shift factor
#' @param pF scaling factor
#' @param pE exponent
#' @param lambda X-ray wavelength in [nm]
#' @return intensity value
#' @examples
#' Yoneda(tthf = 0, af = 1.5, ai = .22,
#'        delta = 7e-5, beta = 5e-7,
#'        surfRms = 1,
#'        pI0 = 0, pF = 1, pE = 3)
Yoneda <- function(tthf, af, ai, delta, beta, surfRms = 0, pI0 = 0, pF = 1, pE = 1, lambda = .154)
{
  ai <- ai / 180 * pi
  af <- af / 180 * pi
  tran.i <- abs(fresnelT(ai, delta, beta, surfRms = surfRms, lambda = .154)) ^ 2
  tran.f <- abs(fresnelT(af, delta, beta, surfRms = surfRms, lambda = .154)) ^ 2
  afT <- alphaT(af, delta, beta)
  aiT <- alphaT(ai, delta, beta)
  Q <- angle2Q(ai = aiT, af = afT, tthf = tthf, lambda = lambda)
  Qabs <- abs(sqrt(Q$Qx ^ 2 + Q$Qy ^ 2 + Q$Qz ^ 2))
  intens <- pI0 + pF * tran.i * tran.f * exp(-pE * Qabs ^ 2)
  return(intens)
}

#' alphaT(alpha, delta, beta)
#' @description calculate angle of refracted beam
#' @author K. Juraic
#' @param ai angle of incidence in [deg.]
#' @param delta Re part of index of refraction
#' @param beta Im part of index of refraction
#' @return at angle of refracted beam in [deg.]
#' @examples
#' alphaT(ai = 0.5, delta = 7e-5, beta = 5e-7)
alphaT <- function(ai, delta, beta){
    at <- sqrt(ai^2 - 2 * delta + 2i * beta)
  return(at)
}

#' Fresnel reflection coeficient for rough interface
#' @description calculate Fresnel coeficient for reflection at rough surface
#' @author K. Juraic
#' @param ai angle of incidence in [rad.]
#' @param delta Re part of index of refraction
#' @param beta Im part of index of refraction
#' @param surfRms RMS surface roughness in [nm]
#' @param lambda X-ray wavelength
#' @return fr reflection coeficient
#' @examples
#' fresnelR(ai = 0.5, delta = 7e-5, beta = 5e-7, surfRms = 0)
fresnelR <- function(ai, delta, beta, surfRms = 0, lambda = .154){
  k0 <- 2 * pi / lambda
  at <- alphaT(ai, delta, beta)
  aiSin <- sin(ai)
  atSin <- sin(at)
  nr <- 1 - 2 * (delta + 1i * beta)
  rmsExp <- exp(-2 * k0^2 * aiSin * atSin * surfRms^2)
  rf <- (aiSin - nr * atSin) / (aiSin + nr * atSin) * rmsExp
  return(rf)
}

#' Fresnel transmission coeficent for rough interface
#' @description calculate Fresnel coeficient for transmission at rough surface
#' @author K. Juraic
#' @param ai angle of incidence in [rad.]
#' @param delta Re part of index of refraction
#' @param beta Im part of index of refraction
#' @param surfRms RMS surface roughness in [nm]
#' @param lambda X-ray wavelength
#' @return ft trnamsission coeficient
#' @examples
#' fresnelT(ai = 0.5, delta = 7e-5, beta = 5e-7, surfRms = 0)
fresnelT <- function(ai, delta, beta, surfRms = 0, lambda = .154){
  k0 <- 2 * pi / lambda
  at <- alphaT(ai, delta, beta)
  aiSin <- sin(ai)
  atSin <- sin(at)
  nr <- 1 - 2 * (delta + 1i * beta)
  rmsExp <- exp(-k0^2 * (aiSin - atSin)^2 * surfRms^2/2)
  tf <- 2 * aiSin / (aiSin + nr * atSin) * rmsExp
  return(tf)
}

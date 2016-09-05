#' penDepth
#' @description penetration depth calculation
#' @author K. Juraic
#' @param ai angle of incidence
#' @param beta Im. part of index of refraction
#' @param delta Re part of index of refraction
#' @param lambda wavelength [nm]
#' @return calculated penetration depth in [nm]
#' @examples
#' ai <- seq(.01,1,by=.001)
#' plot(ai,penDepth(ai, 0.758e-5, 1.755e-7), type='l')
penDepth <- function(ai, delta, beta, lambda = 0.154){
  ai <- ai / 180 * pi # angle in deg.
  n <- 1 - delta + 1i * beta
  ac <- acos(n)
  #ac <- sqrt(2*delta)
  out  <- lambda / (2*sqrt(2)*pi) / (sqrt(sqrt((ai^2-ac^2)^2+4*beta^2)-ai^2+ac^2))
  return(Re(out))
}


#' alphaifCalc
#' @description calculate angle alpha.i + alpha.f from from known gisaxs
#'              parameter calibration (central pixel, sample to detector distance)
#' @author K. Juraic
#' @param pxV pixel position
#' @param pxCenterV  cordinate of direct beam position
#' @param samDet sample to detector distance
#' @param tiltV vertical detector tilt
#' @param pxSizeV detector pixel size in vertical direction
#' @return angle alpha.i + alpha.f for pxV in deg.
#' @examples
#' alphaifCalc(100, 50, 1700, -5.21)
alphaifCalc <- function(pxV, pxCenterV, samDet, tiltV = 0, pxSizeV = .172){
  # Side View( detector left, sample right) => counterclockwise rotation
  pxDist <- (pxV - pxCenterV) * pxSizeV
  #pxL <- sqrt(samDet^2 + pxDist^2 - 2 * samDet * pxDist * cos((90+tiltV)/180*pi))
  #alphaif <- asin(pxDist/pxL * sin(((90+tiltV)/180*pi))) / pi * 180
  alphaif <- atan(pxDist/samDet) / pi * 180
  return(alphaif)
}

#' tthfCalc
#' @description calculate angle 2*theta exit from from known gisaxs
#'              parameter calibration (central pixel, sample to detector distance)
#' @author K. Juraic
#' @param pxH pixel position (horizontal coordinate)
#' @param pxCenterH  cordinate of direct beam position
#' @param samDet sample to detector distance
#' @param tiltH horizontal detector tilt
#' @param pxSizeH detector pixel size in horizontal direction
#' @return angle 2*theta exit for pxH
#' @examples
#' tthfCalc(100, 50, 1700, -5.21)
tthfCalc <- function(pxH, pxCenterH, samDet, tiltH = 0, pxSizeH = .172){
  # view from top -> counter clockwise rotatio
  pxDist <- (pxH - pxCenterH) * pxSizeH
  #pxL <- sqrt(samDet^2 + pxDist^2 - 2 * samDet * pxDist * cos((90+tiltH)/180*pi))
  #tthf <- asin(pxDist/pxL * sin(((90+tiltH)/180*pi))) / pi * 180
  tthf <- atan(pxDist/samDet) / pi * 180
  return(tthf)
}

#' angle2Q
#' @description calculate Q vector from angle pair (alpha.f, tth.f)
#' @author K. Juraic
#' @param ai angle of incidence in spec plane
#' @param af exit angle in spec plane
#' @param tthf 2*theta exit angle
#' @param lambda wavelength
#' @return Q vector(qx,qy,qz)
#' @examples
#' angle2Q(.22, .28, 0)
angle2Q <- function(ai, af, tthf, lambda = .154){
  k0 <- 2 * pi / lambda
  aiRad <- ai / 180 * pi
  afRad <- af / 180 * pi
  tthfRad <- tthf / 180 * pi
  Qx <- k0 * (cos(afRad) * cos(tthfRad) - cos(aiRad))
  Qy <- k0 * (cos(afRad) * sin(tthfRad))
  Qz <- k0 * (sin(aiRad) + sin(afRad))
  return(data.frame(Qx = Qx, Qy = Qy, Qz = Qz))
}

#' Angle of refracted beam calculation
#' @description angle of transmited beam trough surface
#' @author K. Juraic
#' @param ai angle of incidence
#' @param beta Im. index of refracion
#' @param delta Re index of refraction
#' @return alphaT angle of refracted beam
#' @examples
#' alphaTcalc(.22,0.758e-5, 1.755e-7)
alphaTcalc <- function(ai, delta, beta){
  ai <- ai / 180 * pi
  alphaT <- sqrt(ai^2 - 2*delta + 2i * beta)
  alphaT <- alphaT / pi * 180

  return(Re(alphaT))
}



#' Section Angle and Q scale calculation
#' @description calculate Angular and Q scale for section from pixels array
#' @author K. Juraic
#' @param sec section data.frame(pxH, pxV, Imean, Isd)
#' @param ai angle of incidence in deg.
#' @param calib data.frame with calibration parameters (pxcH, pxcV, saDet, tiltH, tiltV)
#' @return data.frame(pxH, pxV, Imean, Isd, aif, tthf, af, Qx, Qy, Qz)
#' @examples
#'    \dontrun{sc <- secCalcAngleQ(secDf, ai, calibDf)}
#'
secCalcAngleQ <- function(sec, ai, calib){
  sec$aif <- alphaifCalc(pxV = sec$pxV, pxCenterV = calib$pxcV, samDet = calib$saDet, tiltV = calib$titlV)
  sec$tthf <- tthfCalc(pxH = sec$pxH, pxCenterH = calib$pxcH, samDet = calib$saDet, tiltH = calib$tiltH)
  sec$af <- sec$aif - ai
  tmpQ <- angle2Q(ai = ai, af = sec$af, tthf = sec$tthf, lambda = 0.154)
  sec$Qx <- tmpQ$Qx
  sec$Qy <- tmpQ$Qy
  sec$Qz <- tmpQ$Qz
  return(sec)
}


# readNikaCalib -----------------------------------------------------------
#' Read Nika calibration parameters
#' @author Krunoslav Juraic
#' @description Readcalibraion parameters calculated by Nika Igor Pro package.
#'              Parameters shoud be written to tab delimited file. First row
#'              shoud have parametr names (pxcH	pxcV	saDet	tiltH	tiltV) and
#'              second row values. If calibration file does not exist program
#'              asks to insert parameters manualy and then write this
#'              parameters to file.
#' @param calibFile full path for calibration file
#' @return data.frame with calibration data
#' @examples
#'        \dontrun{readNikaCalib(calibFile = "calib.txt")}
readNikaCalib <- function(calibFile){
  if (file.exists(calibFile)) {
    calib <- read.table(file = calibFile, header = TRUE, sep = "\t")
  } else {
    cat("File [", calibFile, "does not exist!\n")
    cat("Insert calibration parameters manualy\n")
    calib <- data.frame(pxcH = 0,
                        pxcV = 0,
                        saDet = 0,
                      	tiltH = 0,
	                      tiltV = 0)
    calib <- edit(calib)
    write.table(x = calib, file = calibFile,
                col.names = TRUE, row.names = FALSE,
                sep = "\t", quote = FALSE)
    cat("Calibration parameters written to file:\n\t", calibFile)
  }
  assign(x = "calib", value = calib, envir = .GlobalEnv)
  cat("Calibration data read successfully in calib data.frame\n")
  return(calib)
}


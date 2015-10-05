
# Function: writeP00 ------------------------------------------------------
#' writeP00
#' @description export 1D SAXs data to p00 format
#' @author K. Juraic
#' @param fname file name
#' @param qExp Q vector
#' @param iExp intensity values
#' @param distSD sample to detector distance
#' @param pxC direct beam position in px
#' @param pxSize detector pixel size
#' @param pxN number of pixels
#' @param lambda wavelength
#' @return out data.frame(px,intens)
#' @examples
#' write\dontrun{P00(fname, qExp, iExp, distSD, pxC, pxSize, pxN, lambda)}
#'
writeP00 <- function(fname, qExp, iExp, distSD, pxC, pxSize, pxN, lambda){
  fname <- sub("\\..*",".p00",fname,perl = TRUE)
  px <- 1:pxN
  tth <- atan((px-pxC)*pxSize/distSD)
  q <- 4*pi/lambda * sin(tth/2)
  interpExp <- approx(x = qExp, y = iExp, xout = q, yleft = 0, yright = 0)
  write.table(data.frame(interpExp$y), file = fname, row.names=FALSE, col.names=FALSE)
  cat("Date saved in [",fname,"]\n")
  return(data.frame(px,intens=interpExp$y))
}


# Function: writeIsGISAXSDat ----------------------------------------------
#' writeIsGisaxsDat
#' @description write 1D experimental data in IsGISAXS dat format for fiting
#'              Experimental data should be given in format
#'              data.frame(T, sin(2tth), sin(alphaF), Intensity, sigmaI)
#' @author K. Juraic
#' @param fName IsGISAXS dat file name
#' @param af exit angle alpha.f
#' @param tthf exit angle 2 theta
#' @param intens intensity
#' @param sigmaI intensity error
#' @param fitTF iclud this point i fit (T OR f)
#' @param weight weigth fot evry point of cross section
#' @param scaleFac The scale factor of the cross section with hard bounds below.
#' @param shiftFac The shift factor of the cross section with hard bound below.
#' @param comments experiment or data comment
#' @return fName
#' @examples
#' \dontrun{writeIsGisaxsDat(fName = "exp_data.dat", af, tthf, intens, sigmaI, fitTF,
#' weight = 1, scaleFac = 1, shiftFac = 0, comments = "")}
#'
writeIsGisaxsDat <- function(fName, af, tthf, intens, sigmaI, fitTF,
                             weight = rep(1:length(fName)), scaleFac = 1, shiftFac = 0,
                             comments = ""){

  fName <- sub("\\..*",".dat",fName,perl = TRUE)
  cat(paste(replicate(80, "#"), collapse = ""), file = fName, append = FALSE)
  cat(paste("#", comments), file = fName, append = TRUE)
  cat(paste("#", fName), file = fName, append = TRUE)
  cat("# Weight,   Scale factor,   Shift factor", file = fName, append = TRUE)
  cat(paste(weight, scaleFac, "F", shiftFac, "F"), file = fName, append = TRUE)
  cat(paste("0  1e5          0   1e5"), file = fName, append = TRUE)
  cat(paste("F  F           F  F"), file = fName, append = TRUE)
  cat(paste("# DeltaOmega(deg), DeltaOmega(deg)"), file = fName, append = TRUE)
  cat(paste("0.                0"), file = fName, append = TRUE)
  cat(paste("# Fitted, Sin(2Thetaf), Sin(Alphaf), Intensity, Error bars"), file = fName, append = TRUE)
  ft <- rep("T", times = length(af))
  df <- data.frame(fitTF, stthf = sin(tthf/180*pi), saf = sin(af/180*pi), intens, sigmaI)
  write.table(df, file = fName, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "")
  cat(" Data saved to file [", fName, "]")
}

#' export SAXS data in pdh format
#'
#' @param param_list experiment parameters
#' @param saxs_dat experimental data
#' @param file_name output pdh file name
#'
#' @return none
#' @examples
#' saxs_dat <- data.frame(q = 1:100, I = sqrt(1:100), Ierr = sqrt(1:100)*.1)
#' file_name <- "proba.phd"
#' param_list <- list(title = "AustroSAXS",
#'                    keywords = "SAXS",
#'                    number_of_points = nrow(saxs_dat),
#'                    sample_to_detector_mm = 1741.829,
#'                    normalization_factor = 1,
#'                    lambda_nm = 1.54e-1
#' )
#'    export_pdh(param_list, saxs_dat, file_name)
#'
export_pdh <- function(param_list, saxs_dat, file_name)
{
  cat("Writin to file:", file_name,"\n")
  file_con <- file(file_name)
  pdh3 <- sprintf("%9d %9d %9d %9d %9d %9d %9d %9d",param_list$number_of_points,0,0,0,0,0,0,0)
  pdh4 <- sprintf("%14.6e %14.6e %14.6e %14.6e %14.6e", 0, param_list$sample_to_detector_mm, 0, param_list$normalization_factor, param_list$lambda)
  pdh5 <- sprintf("%14.6e %14.6e %14.6e %14.6e %14.6e", 0, 0, 0, 0, 0)
  phd_header <- c(param_list$title, param_list$keywords, pdh3, pdh4, pdh5)
  saxs_dat_foramted <- sprintf("%14.6e %14.6e %14.6e", saxs_dat[,1],saxs_dat[,2], saxs_dat[,3])
  writeLines(text = c(phd_header, saxs_dat_foramted), con = file_con, sep = '\n', )
  close(file_con)
}



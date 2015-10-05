#' preview horizontal and vertical 1D section for 2D GISAXS image
#' @description at specific pixel and withspecific width preview horizontal and
#'              vertical corssection for 2D GISAXS image
#' @author K. Juraic
#' @param gisaxs_mat matrix with GISAXS pattern
#' @return section_lst last viewed sections
#' @examples
#' \dontrun{section_preview(gisaxs.mat)}
section_preview <- function (gisaxs_mat) {
  px_h <- NULL # avoid NOTE in package check
  px_v <- NULL # avoid NOTE in package check
  px_hwd <- NULL # avoid NOTE in package check
  s <- NULL # avoid NOTE in package check
  mat_dim <- dim(gisaxs_mat)
  manipulate({
      sec_v<-colMeans(gisaxs_mat[(px_v-px_hwd):(px_v+px_hwd),], na.rm = TRUE)
      sec_h<-rowMeans(gisaxs_mat[,(px_h-px_hwd):(px_h+px_hwd)], na.rm = TRUE)
      int_max <- max(c(sec_v,sec_h))
      s_h<-seq(1,mat_dim[1],by=s)
      s_v<-seq(1,mat_dim[2],by=s)
      s_z<- gisaxs_mat[s_h,s_v]
      par(mfrow = c(1,2))
      image(s_h, s_v, log(s_z),col = rainbow(100))
      abline(h=px_h - px_hwd,col='red')
      abline(h=px_h + px_hwd,col='red')
      abline(v=px_v - px_hwd,col='blue')
      abline(v=px_v + px_hwd,col='blue')
      plot(sec_h,type='l',col='red',log='y');
      points(sec_v,type='l',col='blue');
      abline(v=px_h - px_hwd,col='red')
      abline(v=px_h + px_hwd,col='red')
      abline(v=px_v - px_hwd,col='blue')
      abline(v=px_v + px_hwd,col='blue')
    },
    px_h = slider(min = 1, max = mat_dim[1], initial = mat_dim[1]/2),
    px_v = slider(min = 1, max = mat_dim[2], initial = mat_dim[2]/2),
    px_hwd = slider(min = 1, max = 50, initial = 1),
    s = picker("1x1" = 1,"2x2" = 2,"4x4" = 4,"8x8" = 8, "16x16" = 16, initial = "16x16")
  )
}


#' sectionV(mat, pos, width = 10)
#' @description make a vertical section of 2D matrix at specific position and width
#' @author K. Juraic
#' @param mat 2D matrix
#' @param pos  trace position
#' @param width trace width
#' @return out.df section data.frame(pxH, pxV, Imean, Isd)
#' @examples
#' \dontrun{sectionH(gisaxs.mat)}
sectionV <- function(mat, pos, width = 10){
  pxH <- rep(pos,times = dim(mat)[2])
  pxV <- 1:dim(mat)[2]
  Imean <- colMeans(mat[(pos - width):(pos + width),], na.rm = TRUE)
  Isd <- apply(mat[(pos - width):(pos + width),], 2, sd, na.rm = TRUE)
  out.df <- data.frame(pxH, pxV, Imean, Isd)
  return(out.df)
}

#' sectionH(mat, pos, width = 10)
#' @description make a horizontal section of 2D matrix at specific position and width
#' @author K. Juraic
#' @param mat 2D matrix
#' @param pos  trace position
#' @param width trace width
#' @return out.df section data.frame(pxH, pxV, Iraw)
#' @examples
#' \dontrun{sectionH(gisaxs.mat)}
sectionH <- function(mat, pos, width = 10){
  pxH <- 1:dim(mat)[1]
  pxV <- rep(pos,times = dim(mat)[1])
  Imean <- rowMeans(mat[,(pos - width):(pos + width)], na.rm = TRUE)
  Isd <- apply(mat[,(pos - width):(pos + width)], 1, sd, na.rm = TRUE)
  out.df <- data.frame(pxH, pxV, Imean, Isd)
  return(out.df)
}

#' Vertical section for list of tif images
#' @description make a vertical section for list of 2D matrix at specific
#'              position and width. And save to file
#' @author K. Juraic
#' @param tifLst list of 2D matrices
#' @param tifName names of tif files
#' @param pos  trace position in pixels
#' @param width trace width in pixels
#' @param write2file store sections to files
#' @param path folder where to store sections
#' @return sec list of section list(data.frame(pxH, pxV, Imean, Isd))
#' @examples
#' \dontrun{sectionVlst(tifLst, tifName, pos, width)}
sectionVlst <- function(tifLst, tifName, pos, width = 10, write2file = FALSE, path = "./"){
  sec <- list()
  for (i in 1:length(tifLst)) {
    secName <- paste("secV_", pos, "_", width, sub(".tif", ".dat", tifName[i]), sep = "")
    cat("[", i,"/", length(tifLst),"]", secName, "\n")
    sec[[secName]] <- sectionV(mat = tifLst[[i]], pos = pos, width = width)
    if (write2file) {
      write.table(sec[[i]], file = paste(path,secName,sep=""),
                  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  return(sec)
}

#' Horizontal section for list of tif images
#' @description make a horizontal section for list of 2D matrix at specific
#'              positions and width. And save to file
#' @author K. Juraic
#' @param tifLst list of 2D matrices
#' @param tifName names of tif files
#' @param pos  array of trace position in pixels
#' @param width trace width in pixels
#' @param write2file store sections to files
#' @param path folder where to store sections
#' @return sec list of section list(data.frame(pxH, pxV, Imean, Isd))
#' @examples
#' \dontrun{sectionHlst(tifLst, tifName, pos, width)}
sectionHlst <- function(tifLst, tifName, pos, width = 10, write2file = FALSE, path = "./"){
  sec <- list()
  pos <- round(pos)
  for (i in 1:length(tifLst)) {
    secName <- paste("secH_", pos[i], "_", width, sub(".tif", ".dat", tifName[i]), sep = "")
    cat("[", i,"/", length(tifLst),"]", secName, "\n")
    sec[[secName]] <- sectionH(mat = tifLst[[i]], pos = pos[i], width = width)
    if (write2file) {
      write.table(sec[[i]], file = paste(path, secName, sep = ""),
                  col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
    }
  }
  return(sec)
}


#' secLstShow(secLst, titles, ...)
#' @description Show list of matrix sections
#' @author K. Juraic
#' @param secLst list of matrix sections
#' @param titles main title for plots
#' @param colSelect  columns to plot
#' @param ... other parameters log='y'
#' @examples
#' \dontrun{secLstShow(dat, fNames, log='y')}
secLstShow <- function(secLst, titles, colSelect = c(1,3),...){
  id <- NULL # avoid NOTE in package check
  manipulate(
    plot(secLst[[id]][,colSelect[1]],secLst[[id]][,colSelect[2]], main = titles[id], type = 'o', pch = 16, cex = .5, col = 'red', ...),
    id = slider(min = 1, max = length(secLst))
  )
}


#' Add angular scale to section
#'
#' @description Add angular scale columns (aif, af, tthf) to section data frame
#' @author K. Juraic
#'
#' @param sec section data frame (pxH, pxV, Imean, Isd)
#' @param ai angle of incidence in deg.
#' @param calib calibration parameters (pxcH, pxcV, saDet, tiltH, tiltV)
#'
#' @return data.frame(pxH, pxV, Imean, Isd, aif, ai, tthf)
#'
#' @examples
#'  \dontrun{secDfAddAngle(sec, ai, calib)}
secDfAddAngle <- function(sec,ai,calib){
  sec$aif <- alphaifCalc(pxV = sec$pxV, pxCenterV = calib$pxcV, samDet = calib$saDet, tiltV = calib$tiltV)
  sec$tthf <- tthfCalc(pxH = sec$pxH, pxCenterH = calib$pxcH, samDet = calib$saDet, tiltH = calib$tiltH )
  sec$af <- sec$aif - ai
  return(sec)
}


#' Add Q scale to section
#'
#' @description Add Q scale columns (qx, Qy, Qz, Qabs) to section data frame
#' @author K. Juraic
#' @param sec section data frame (pxH, pxV, Imean, Isd, aif, ai, tthf)
#' @param ai angle of incidence in deg.
#' @param lambda wavelength in nm
#'
#' @return data.frame(pxH, pxV, Imean, Isd, aif, ai, tthf, Qx, Qy, Qz, Qabs)
#'
#' @examples
#'  \dontrun{secDfAddAngle(sec, ai, lambda)}
secDfAddQ <- function(sec, ai, lambda = .154){
  tmpQ <- angle2Q(ai = ai, af = sec$af, tthf = sec$tthf, lambda = lambda)
  sec$Qx <- tmpQ$Qx
  sec$Qy <- tmpQ$Qy
  sec$Qz <- tmpQ$Qz
  sec$Qabs <- sqrt(sec$Qx ^ 2 + sec$Qy ^ 2 + sec$Qz ^ 2)
  if (sec$Qy < 0 | sec$Qz < 0)
    sec$Qabs <- -sec$Qabs
  return(sec)
}


#' Add angular scale to list of section
#'
#' @description Add angular scale columns (aif, af, tthf) to list of section data frame
#' @author K. Juraic
#'
#' @param secLst list of section data frames (pxH, pxV, Imean, Isd)
#' @param aiLst list of angle of incidence in deg.
#' @param calib calibration parameters (pxcH, pxcV, saDet, tiltH, tiltV)
#'
#' @return seclst list of section data.frame(pxH, pxV, Imean, Isd, aif, ai, tthf)
#'
#' @examples
#'  \dontrun{lstSecDfAddAngle(secLst, aiLst, calib)}
lstSecDfAddAngle <- function(secLst, aiLst, calib){
  for (id in 1:length(secLst)) {
    secLst[[id]] <- secDfAddAngle(secLst[[id]],aiLst[id],calib = calib)
  }
  return(secLst)
}


#' Add Q scale to list of section
#'
#' @description Add Q scale columns (qx, Qy, Qz) to list of section data frame
#' @author K. Juraic
#' @param secLst list of section data frame (pxH, pxV, Imean, Isd, aif, ai, tthf)
#' @param aiLst list of angle of incidence in deg.
#' @param lambda wavelength in nm
#'
#' @return secLst list of data.frame(pxH, pxV, Imean, Isd, aif, ai, tthf, Qx, Qy, Qz)
#'
#' @examples
#'  \dontrun{lstSecDfAddAngle(secLst, aiLst, lambda)}
lstSecDfAddQ <- function(secLst, aiLst, lambda = .154){
  for (id in 1:length(secLst)) {
    secLst[[id]] <- secDfAddQ(secLst[[id]],aiLst[id], lambda = lambda)
  }
  return(secLst)
}


#' save section data.frame to file
#'
#' @param sec section data frame
#' @param fname section file name (secH_430_10_sample_name.dat)
#' @param colSubset subset of column to save to file

#' @examples
#'  \dontrun{sec2file(sec, fname)}
sec2file <- function(sec, fname, colSubset = 1:nrow(sec)){
  sec <- sec[!is.na(sec$Imean),] # remove NA intensity points
  sec <- sec[,colSubset]
  write.table(x = sec,
              file = fname,
              sep = "\t",
              quote = FALSE,
              row.names = FALSE,
              col.names = TRUE)
  cat("Section writen to file:", fname, "\n")
}


#' save list of section data.frame to file
#'
#' @param secLst list of section data frame
#' @param colSubset subset of column to save to file

#' @examples
#'  \dontrun{lstSec2file(sec, fname)}
lstSec2file <- function(secLst, colSubset = 1:nrow(secLst[[1]])){
  for (id in 1:length(secLst))
    sec2file(sec = secLst[[id]], fname = names(secLst)[id], colSubset = colSubset)
}


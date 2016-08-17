
#' position of maximum in part of matrix
#'
#' @param mat 2D matrix
#' @param subset_H horizontal subsection
#' @param subset_V vertical subsection
#' @param show display m
#'
#' @return c(h_max, v_max) cordinates of matrix
#'
#' @examples
#'  \dontrun{mat_max_pos(mat)}
mat_max_pos <- function(mat, subset_H = 1:dim(mat)[1], subset_V = 1:dim(mat)[2], show = TRUE){
  #require(rpeakfit)
  dim_mat <- dim(mat)
  mat_subset <- mat[subset_H,subset_V]
  image(subset_H,subset_V,log(mat_subset), col = rainbow(100))
  max_pos <- which(mat_subset == max(mat_subset, na.rm = TRUE), arr.ind = TRUE)
  scH <- sectionH(mat = mat_subset, pos = max_pos[2], width = 10)
  scV <- sectionV(mat = mat_subset, pos = max_pos[1], width = 10)
  max_pos <- c(subset_H[max_pos[1]],subset_V[max_pos[2]])
  max_pos
  abline(v = max_pos[1], col = 1)
  abline(h = max_pos[2], col = 1)
  cat("Click on image to continue ...\n")
  tmp <- locator(1)
  plot(subset_H,scH$Imean,type = 'l', col = 'red')
  fitH <- fit_gaussianAmpl(x = subset_H, scH$Imean)
  locator(1)
  plot(subset_V,scV$Imean,type = 'l',col = 'blue')
  fitV <- fit_gaussianAmpl(x = subset_V, scV$Imean)
  cat("Click on image to continue ...\n")
  tmp <- locator(1)
  image(1:dim_mat[1], 1:dim_mat[2],log(mat), col = rainbow(100))
  abline(v = max_pos[1], col = 1)
  abline(h = max_pos[2], col = 1)
  max_pos[1] <- coef(fitH)['xc']
  max_pos[2] <- coef(fitV)['xc']
  return(max_pos)
}


#' matrix crop
#' @description crop matrix to specific region
#' @author K. Juraic
#' @param mat matrix
#' @param bndRow c(px1,px2) horizontal boundaries for croping in pixels
#' @param bndCol c(px1,px2) vertical boundaries for croping in pixels
#' @return list(pxRow, pxCol, mat) cropped matrix with croped horizontal and vertical scale
#' @examples
#' m <- matrix(1:100,nrow = 10, ncol = 10, byrow = TRUE)
#' crop_px_h <- c(2,8)
#' crop_px_v <- c(3,6)
#' matCrop(m,crop_px_h,crop_px_v)
matCrop <- function(mat,bndRow,bndCol)
{
  cropRow <- bndRow[1]:bndRow[2]
  cropCol <- bndCol[1]:bndCol[2]
  mat <- mat[cropRow,cropCol]
  matLst <- list(pxRow=cropRow,pxCol=cropCol,mat=mat)
  return(matLst)
}

  #' matrix shift
  #' @description shift matrix rows and columns
  #' @author K. Juraic
  #' @param m matrix
  #' @param shift c(h,v) number of columns and rows for shift
  #' @return m_shift shifted matrix
  #' @examples
  #' m <- matrix(1:100,nrow = 10, ncol = 10, byrow = TRUE)
  #' shift <- c(0,-2)
  #' matrix_shift(m,shift)
matrix_shift <- function (m, shift) {
  # horizontal shift
  m_shift <- m
  if (shift[2]>0) {
    m_cut <- m_shift[,1:(dim(m)[2]-shift[2])]
    m_add <- matrix(m_shift[,1],nrow = dim(m)[1], ncol = shift[2], byrow = FALSE)
    m_shift <- cbind(m_add,m_cut)
  }
  if (shift[2]<0) {
    m_cut <- m_shift[,abs(shift[2]-1):dim(m)[2]]
    m_add <- matrix(m_shift[,dim(m)[2]],nrow = dim(m)[1], ncol = abs(shift[2]), byrow = FALSE)
    m_shift <- cbind(m_cut,m_add)
  }
  # vertical shift
  if(shift[1]>0){
    m_cut <- m_shift[(shift[1]+1):(dim(m)[1]),]
    m_add <- matrix(m_shift[dim(m)[1],],nrow = shift[1], ncol = dim(m)[2], byrow = TRUE)
    m_shift <- rbind(m_cut,m_add)
  }
  if (shift[1]<0) {
    m_cut <- m_shift[1:(dim(m)[1]-abs(shift[1])),]
    m_add <- matrix(m_shift[1,],nrow = abs(shift[1]), ncol = dim(m)[2], byrow = TRUE)
    m_shift <- rbind(m_add,m_cut)
  }

  m_shift
}


#' matLstShow(matLst, titles, ...)
#' @description Show list of image matrix data
#' @author K. Juraic
#' @param matLst list of matrix data
#' @param titles main title for plots
#' @param ... other parameters log='z'
#' @examples
#' \dontrun{matLstShow(dat, fNames, log='z')}
matLstShow <- function(matLst, titles, ...){
  id <- NULL # ingnore NOTE in package check
  sampling <- NULL # ignore NOTE in package check
  manipulate(
    matPlot(matLst[[id]], main = titles[id], ...),
    id = slider(min = 1, max = length(matLst)),
    sampling = slider(min = 1, max = 32, initial = 32)
  )
}



#' gisaxs2D_angle(mat, aif, tthf, tif_name, save2file = FALSE, width = 600, heiht = 400, title = TRUE)
#' @description display 2D matrix data with Q scale and/or save in png file
#' @author K. Juraic
#' @param mat z values (intensity)
#' @param aif alpha_if scale values
#' @param tthf tthf scale values
#' @param tif_name tif file name
#' @param save2file save to file
#' @param width image width
#' @param height image height
#' @param title to put title or not
#' @examples
#' \dontrun{gisaxs2D_Q(mat, aif, tthf, tif_name,
#'                      save2file = FALSE,
#'                      width = 600, heiht = 400, title = TRUE)}
gisaxs2D_angle <- function(mat, aif, tthf, tif_name, save2file = FALSE, width = 2, height = 1, title = TRUE){
  xlab <- expression(bold(paste(2*theta[f]," [deg.]")))
  ylab <- expression(bold(paste(alpha[i] + alpha[f]," [deg.]")))
  clab <- 'Intensity [a.u.]'
  x <- tthf
  y <- aif
  #z <- log(mat)
  z <- mat
  main_title <- ""
  if(title)
    main_title <- tif_name
  image2D(x = x,
          y = y,
          z = z,
          asp = 1,
          NAcol = "black",
          col = palette(rich.colors(255)),
          #zlim = c(1,max(as.vector(z),na.rm = TRUE)),
          zlim = c(median(mat, na.rm = TRUE)*.9,max(as.vector(mat),na.rm=TRUE)),
          log = "z",
          xlab = xlab,
          ylab = ylab,
          clab = clab,
          main = main_title,
          font.lab = 2,
          font.axis = 2,
          cex.lab = 1.6,
          mgp = c(1.8,.5,0)
  )
  if (save2file) {
    (png_name <- sub(".tif", "_angle.png", tif_name))
    png(filename = png_name, width = width, height = height)
    image2D(x = x,
            y = y,
            z = z,
            asp = 1,
            NAcol = "black",
            col = palette(rich.colors(255)),
            zlim = c(median(mat, na.rm = TRUE)*.9,max(as.vector(mat),na.rm=TRUE)),
            log = 'z',
            xlab = xlab,
            ylab = ylab,
            clab = clab,
            main = main_title,
            font.lab = 2,
            font.axis = 2,
            cex.lab = 1.6,
            mgp = c(1.8,.5,0)
    )
    dev.off()
  }
}

#' gisaxs2D_Q(mat, Qy, Qz, tif_name, save2file = FALSE, width = 600, heiht = 400, title = TRUE)
#' @description display 2D matrix data with Q scale and/or save in png file
#' @author K. Juraic
#' @param mat z values (intensity)
#' @param Qy Qy scale values
#' @param Qz Qz scale values
#' @param tif_name tif file name
#' @param save2file save to file
#' @param width image width
#' @param height image height
#' @param title to put title or not
#' @examples
#' \dontrun{gisaxs2D_Q(mat, Qy, Qz, tif_name, save2file = FALSE,
#'                      width = 600, heiht = 400, title = TRUE)}
gisaxs2D_Q <- function(mat, Qy, Qz, tif_name, save2file = FALSE, width = 600, height = 400, title = TRUE){
  xlab <- expression(bold(paste(Q[y]," [",nm^{-1},"]")))
  ylab <- expression(bold(paste(Q[z]," [",nm^{-1},"]")))
  clab <- 'Intensity [a.u.]'
  x <- Qy
  y <- Qz
  z <- log(mat)
  main_title <- ""
  if(title)
    main_title <- tif_name
  op <- par()
  par()
  image2D(x = Qy,
          y = Qz,
          z = z,
          asp = 1,
          NAcol = "black",
          col = palette(rich.colors(255)),
          zlim = c(1,max(as.vector(z),na.rm = TRUE)),
          xlab = xlab,
          ylab = ylab,
          clab = clab,
          main = main_title,
          font.lab = 2,
          font.axis = 2,
          cex.lab = 1.6,
          mgp = c(1.8,.5,0)
  )
  if (save2file) {
    (png_name <- sub(".tif", "_Q.png", tif_name))
    png(filename = png_name, width = width, height = height)
    image2D(x = x,
            y = y,
            z = z,
            asp = 1,
            NAcol = "black",
            col = palette(rich.colors(255)),
            zlim = c(1,max(as.vector(z),na.rm = TRUE)),
            xlab = xlab,
            ylab = ylab,
            clab = clab,
            main = main_title,
            font.lab = 2,
            font.axis = 2,
            cex.lab = 1.6,
            mgp = c(1.8,.5,0)
    )
    dev.off()
  }
}


#' matPlot(z, x = 1:dim(z)[1], y = 1:dim(z)[2], scaleType = '', sampling = 1, ...)
#' @description display 2D matrix data with specific scale
#' @author K. Juraic
#' @param z z values (intensity)
#' @param x x scale values
#' @param y y scale values
#' @param scaleType 'q'= Q scale, 'a'=angular scale
#' @param sampling skip sampling pixels for fast drawing
#' @param ... other parameters log='z', main='title'
#' @examples
#' m <- matrix(1:100,nrow = 10, ncol = 10, byrow = TRUE)
#' matPlot(m)
matPlot <- function(z, x = 1:dim(z)[1], y = 1:dim(z)[2], scaleType = '', sampling = 1, ...){
  xind<-seq(1,length(x),by=sampling)
  yind<-seq(1,length(y),by=sampling)
  samplex<-x[xind]
  sampley<-y[yind]
  samplez<-z[xind,yind]
  if (scaleType == 'q') {
    # Q scale
    xlab <- expression(bold(paste(Q[H]," [",nm^{-1},"]")))
    ylab <- expression(bold(paste(Q[V]," [",nm^{-1},"]")))
    clab <- 'Intensity [a.u.]'
  } else if (scaleType == 'a') {
    # Angular scale
    xlab <- expression(bold(paste(2*theta[f]," [",nm^{-1},"]")))
    ylab <- expression(bold(paste(alpha[f]," [",nm^{-1},"]")))
    clab <- 'Intensity [a.u.]'
  } else {
    xlab <- ''
    ylab <- ''
    clab <- ''
  }
  image2D(x = samplex,
          y = sampley,
          z = samplez,
          asp = 1,
          NAcol = "black",
          col=palette(rich.colors(255)),
          zlim = c(1,max(as.vector(z),na.rm=TRUE)),
          xlab = xlab,
          ylab = ylab,
          clab = clab,
          font.lab=2,
          font.axis=2,
          ...
  )
}

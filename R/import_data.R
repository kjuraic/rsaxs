
#' read list of fit2D chi file
#' @description read fit2D integrated SAXS data
#' @author K. Juraic
#' @param lstFname list of file names
#' @return data.frame(q,intens)
#' @examples
#' \dontrun{chi_lst <- c("1.chi","2.chi")
#' readLstFit2dChi(chi_lst)
#' }
readLstFit2dChi <- function(lstFname)
{
  dat <- list()
  for(i in 1:length(lstFname))
    dat[[i]] <- readFit2dChi(lstFname)
  return(dat)
}

#' readNika1d
#' @description read integrated SAXS 1D data reduced by NIKA IgorPro package
#' @author K. Juraic
#' @param fname data file name
#' @return data.frame(q, intens, err, unknown)
#' @examples
#' \dontrun{readNika1d("saxs_data.dat")}
#'
readNika1d <- function(fname)
{
  dat <- read.table(fname,skip=23,col.names=c("q","intens","err","unknown"))
  return(dat)
}

#' Read Nika section file
#' @author K. Juraic
#' @description read horizontal/vertical section od 2D SAXS images done by Igor Pro Nika package
#' @param fname Nika section file name
#' @return data.frame(Q, Qy, Qz, Intensity, Error)
#' @examples
#'    \dontrun{read_nika_sec("nika_sec_HLp_0.dat")}
read_nika_sec <- function(fname){
  tmp <- read.table(fname,skip = 27, col.names = c("Q", "Qy", "Qz", "Intensity", "Error"))
  return(tmp)
}


#' readLstNika1d
#' @description read list of NIKA IgorPro integrated SAXS data
#' @author K. Juraic
#' @param lstFname list of file names
#' @return data.frame(q,intens)
#' @examples
#' \dontrun{saxs_chi_lst <- c("1.chi","2.chi")
#' readLstNika1d(saxs_chi_lst)
#' }
readLstNika1d <- function(lstFname)
{
  dat <- list()
  for(i in 1:length(lstFname))
    dat[[i]] <- readNika1d(lstFname)
  return(dat)
}


#' readNika1d_waxs
#' @description read integrated WAXS 1D data reduced by NIKA IgorPro package
#' @author K. Juraic
#' @param fname data file name
#' @return data.frame(tth, Imean, Ierr, unknown)
#' @examples
#' \dontrun{readNika1d_waxs("waxs_data.dat")}
#'
readNika1d_waxs <- function(fname)
{
  dat <- data.frame(tth = NA, Imean = NA, Ierr = NA, unknown = NA)
  cat("Reading file:",fname," ...\n")
  if (file.exists(fname)) {
    dat <- read.table(fname,skip = 23,col.names = c("tth","Imean","Ierr","unknown"))
    cat("\t[OK!]\n")
  } else {
    cat("\t[File not exist!]\n")
  }
  return(dat)
}


# Function readChi --------------------------------------------------------
#' readFit2dChi
#' @description read chi file (fit2D integrated SAXS data)
#' @author K. Juraic
#' @param chiName chi file name
#' @return datOut data.frame(Q[nm], Intensity[a.u])
#' @examples
#'  \dontrun{readFit2dChi("saxs_data.chi")}
#'
readFit2dChi <- function(chiName = file.choose()){
  if(file.exists(chiName)){
    datOut <- read.table(file = chiName, skip = 4, col.names = c("qpernm", "intensity"))
  } else {
    cat("File does not exist!!!\n")
  }
  return(datOut)
}


# Function readLogBook ----------------------------------------------------
#' readLogBook( logFile)
#' @description read AustroSAXS logbook file. Call without parameter will open file.chose dialog
#' @author K. Juraic
#' @param logFile logbook file name,
#' @return out data.frame(name, angle, time, comment)
#' @examples
#'  \dontrun{readLogBook("saxs_data.log")}
#'
readLogBook <- function (logFile = file.choose()) {
  smpLog <- read.table(file = logFile, header = TRUE, sep = '\t')
  smpLog[,1] <- as.character(smpLog[,1])
  #smpLog$comment <- as.character(smpLog$comment)
  cat('Number of samples = ', nrow(smpLog))
  return(smpLog)
}

# Function writeLogBook ----------------------------------------------------
#' writeLogBook( logFile)
#' @description write AustroSAXS logbook file.
#' @author K. Juraic
#' @param logDf logbook data.frame
#' @param logFile logbook file name,
#' @return logbook file name
#' @examples
#'  \dontrun{writeLogBook(logDf)}
writeLogBook <- function(logDf, logFile = file.choose()) {
  write.table(logDf, file = logFile,
              sep = '\t',
              row.names = FALSE, col.names = TRUE,
              quote = FALSE)
  return(logFile)
}

# constructLogBook --------------------------------------------------------
#' Read and insert data in GISAXS/GIWAXS logBook file
#' @author Krunoslav Juraic
#' @description Read existing logBook files and insert additional data if
#'              neceseary. If file does not exist read all tif files from
#'              dataDir and manualy insert aditionay data. Refreshed logBook
#'              data.frame will be saved to file logFile. Function can
#'              automaticaly read exposure times from tif files. Data.frame
#'              logDf is dynamicaly added to global environment
#' @param dataDir folder with tif files
#' @param logFile name of logBook file
#' @return logBook data.frame(name, time, angle, sample)
#' @examples
#' \dontrun{constructLogBook(dataDir = data_path_raw, logFile = logFile)}
constructLogBook <- function(dataDir, logFile = file.choose()){
  setwd(dataDir)
  if (file.exists(logFile)){
    logDf <- readLogBook(logFile = logFile)
    if (!is.element("sample", names(logDf)))
      logDf$sample = rep("", times = nrow(logDf))
    if (!is.element("time", names(logDf)))
      logDf$sample = rep(0, times = nrow(logDf))
  } else {
    fnms <- as.character(list.files(path = dataDir, pattern = "*.tif"))
    logDf <- data.frame(name  = fnms,
                        time  = rep(0, times = length(fnms)),
                        angle = rep(0, times = length(fnms)),
                        sample = rep("", times = length(fnms)))

    expTime <- plyr::aaply(as.character(logDf$name),
                     .margins = 1,
                     .fun = read_pilatus_expTime)
    logDf$time <- expTime
  }
  logDf <- edit(logDf)
  writeLogBook(logDf = logDf, logFile = logFile)
  assign(x = "logDf", value = logDf, envir = .GlobalEnv)
  cat("Logbook write successfull in dataframe logDf!\nUpdated file",
      logFile, "\n")
  return(logDf)
}



# Function: readPilatus1M -------------------------------------------------
#' readPilatus1M
#' @description read Pilatus 1M tif file
#' @author K. Juraic
#' @param fname tiff file name
#' @return mat image data in 2D matrix
#' @examples
#' \dontrun{readPilatus1M("saxs_data.tif")}
#'
readPilatus1M<-function(fname = file.choose())
{
  zz <- file(fname, "rb")
  smece<-readBin(zz, raw(),4096) # header at the begining of tif file
  slika<-readBin(zz,what=integer(),n=981*1043,size=4)
  close(zz)
  m<-matrix(as.integer(slika),nrow=981,ncol=1043,byrow=FALSE)
  m<-m[,1043:1]
  m[m < 0]  <- NA
  return(m)
}

# Function: readPilatus1Mlst -------------------------------------------------
#' readPilatus1Mlst
#' @description read list of Pilatus 1M tif file
#' @author K. Juraic
#' @param fnames tiff file name list
#' @return list(mat) list of image data in 2D matrix
#' @examples
#' \dontrun{readPilatus1Mlist(list)}
#'
readPilatus1Mlst<-function(fnames)
{
  n.images<-length(fnames)
  lst.mat<-list()
  cat(sprintf("Number of images = %d\n",n.images))
  for(i in 1:n.images){
    cat(sprintf("\t[%d/%d] %s\n",i,n.images,fnames[i]))
    lst.mat[[fnames[i]]]<-readPilatus1M(fnames[i])
  }
  return(lst.mat)
}


# Function: read_pilatus_expTime-------------------------------------------------
#' read_pilatu_expTime
#' @description read exposure time directly from Pilatus tif image
#' @author K. Juraic
#' @param tifName tiff file name
#' @return exposure time in seconds
#' @examples
#' \dontrun{read_pilat_expTimeu("saxs_data.tif")}
#'
read_pilatus_expTime <- function(tifName){
  expTime <- NA
  if (file.exists(tifName)) {
    con <- file(tifName)
    tmp <- readLines(con = con, n = 18, skipNul = TRUE)
    tmp <- tmp[-c(1:3)] # Firsr 3 rows not important
    expTime <- as.numeric(strsplit(x = tmp[grepl("Exposure_time",tmp)], split = " ")[[1]][3])
    close(con)
    cat("[", tifName, "] Exp. Time =", expTime, "s\n")
  } else {
    cat("File [", tifName, "] does'n exist!\n")
  }
  return(expTime)
}

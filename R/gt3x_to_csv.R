#' gt3x_to_csv
#'
#' @description 'gt3x_to_csv' converts old ActiGraph GT3X files to CSV format
#' 
#' This script converts ActiGraph GT3X accelerometer files to CSV format using the read.gt3x package.
#' 
#' @details The script automates the conversion process for all GT3X files in a specified directory. After converting each GT3X file to CSV format, it optionally compresses the CSV files using gzip.   # This code was kindly provided by Hansjorg Baurecht (Regensburg).
#' 
#' @param path The directory path where GT3X files are located.
#' @param gzip Logical indicating whether to gzip the resulting CSV files (default = TRUE).
#' 
#' @import DescTools
#' @import data.table
#' @import Rfast
#' @import remotes
#' @import read.gt3x
#' 
#' @export

#===================================
gzip = FALSE #as.logical(commandArgs(TRUE)[2])

gt3x_to_csv <- function(path, gzip = T){
  gt3x <- read.gt3x(path = path, imputeZeroes = TRUE, asDataFrame = FALSE, verbose = TRUE)
  outpath<-gsub(".gt3x", ".csv", path, fixed = T)
  # Convert date formats
  start.time <- unlist(strsplit(as.character(attributes(gt3x)$header$'Start Date'), split = " "))
  start.date <- format(attributes(gt3x)$header$'Start Date', format = "%d.%m.%Y")
  download.time <- unlist(strsplit(as.character(attributes(gt3x)$header$'Download Date'), split = " "))
  download.date <- format(attributes(gt3x)$header$'Download Date', format = "%d.%m.%Y")
  # Write header to csv-file
  header <- paste("------------ Data File Created By ActiGraph ",attributes(gt3x)$header$'Device Type'," R read.gt3x v", packageVersion("read.gt3x"), " Firmware v", attributes(gt3x)$header$'Firmware', " date format dd.MM.yyyy at ", attributes(gt3x)$header$'Sample Rate', " Hz Filter Normal -----------",sep="")
  header <- c(header,paste("Serial Number:",attributes(gt3x)$header$'Serial Number',sep=" "))
  header <- c(header,paste("Start Time",start.time[2]))
  header <- c(header,paste("Start Date",start.date))
  epoch <- ifelse(is.null(attributes(gt3x)$header$'Epoch Period'),"",attributes(gt3x)$header$'Epoch Period')
  header <- c(header,paste("Epoch Period (hh:mm:ss)",epoch))
  header <- c(header,paste("Download Time",download.time[2]))
  header <- c(header,paste("Download Date",download.date))
  addr <- ifelse(is.null(attributes(gt3x)$header$'Current Memory Address'),"",attributes(gt3x)$header$'Current Memory Address')
  header <- c(header,paste("Current Memory Address:",addr))
  voltage <- gsub(",",".",attributes(gt3x)$header$'Battery Voltage')
  mode <- ifelse(is.null(attributes(gt3x)$header$'Mode'),"",attributes(gt3x)$header$'Mode')
  header <- c(header,paste("Current Battery Voltage: ",voltage," Mode=",mode,sep=""))
  header <- c(header,"--------------------------------------------------")
  fileConn <- file(outpath)
  writeLines(header, fileConn)
  close(fileConn)
  # Write gravity data
  outdat <- as.matrix(gt3x)
  # Identify lines with 0, 0, 0 and replace them by NA
  correct1 <- rowsums(outdat)
  correct2 <- rowMins(outdat, value = T)
  outdat[, 1]<- ifelse(correct1 == 0 & correct2 == 0, NA, outdat[, 1])
  outdat[, 2]<- ifelse(correct1 == 0 & correct2 == 0, NA, outdat[, 2])
  outdat[, 3]<- ifelse(correct1 == 0 & correct2 == 0, NA, outdat[, 3])
  # If first line begins with vector of zeros, leaf it as vector of zeros
  if(correct1[1] == 0 & correct2[1] == 0) {outdat[1, ] <- c(0, 0, 0)}
  # replace NA by last value carried forward
  outdat[, 1] <- LOCF(outdat[, 1])
  outdat[, 2] <- LOCF(outdat[, 2])
  outdat[, 3] <- LOCF(outdat[, 3])
  print(paste0("outpath ", outpath))
  fwrite(x = outdat, file = outpath, col.name = F, append = T)
  if(gzip == T){
    command <- paste("gzip", outpath)
    system(command)
  }
}

filenames = dir(path = path, full.names = T)
filenames = grep(x = filenames,pattern = "gt3x", value = T)

for (i in 1:length(filenames)) {
  print(filenames[i])
  gt3x_to_csv(filenames[i])
}

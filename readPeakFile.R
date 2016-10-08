##' @title readPeakFile
##' @description Read bed files from file directory as GRanges objects
##' @return GRanges objects

readPeakFiles <- function(peakFolder, verbose=FALSE, ...) {
  # input param checking
  if(missing(peakFolder)) {
    stop("input param is missing!")
  }
  stopifnot(length(peakFolder)>=1)
  files <- list.files(peakFolder, full.names = TRUE, "\\.bed$")
  f.read <- setNames(
    lapply(files, function(ele_) {
      out <- as(import.bed(ele_), "GRanges")
    }), tools::file_path_sans_ext(basename(files))
  )
  res <- f.read
  return(res)
}

##' @examples myData <- readPeakFiles(peakFolder = "data/")

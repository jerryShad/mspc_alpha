##' @title 
##' @description 
##' @return 
##' @example

.denoise_peakFiles <- function(peakFolder, tau.w=1.0E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(tau.w))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  peaks_rmnoi <- lapply(peakFolder, function(ele_) {
    if(is.null(ele_$p.value)) {
      ele_ <- .pvalueConversion(ele_, pvalueBase = 1L)
    } else {
      res <- subset(ele_, ele_$p.value <= tau.w)
    }
  })
  return(peaks_rmnoi)
}

##' @example total.ERs <- .denoise_peakFiles(peakFolder = myData, tau.w = 1.0E-04)

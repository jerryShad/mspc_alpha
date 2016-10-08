##' @title filtByHitIdx
##' @description 
##' @return 

.filtByHitIdx <- function(peakset, ovHit, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  ovNum <- Reduce('+', lapply(ovHit, lengths))
  keepMe <- sapply(ovNum, function(elm) {
    res <- elm >= min.c
    res
  })
  return(keepMe)
}

#' @example
keep_me <- .filtByIndx(peakset = total.ERs, .hit_3,
                       replicate.type = "Biological")

library(purrr)

#' @example 
keepList <- map(.hit_1, ~.[keep_me])
dropList <- map(.hit_1, ~.[!keep_me])

##' @description peaks that did not meet the sufficient overlap condition won't be proceed to next,
##' instead we keep them as GRanges to get overall discarded peaks for statistical result at the end

Discard.peaks_0 <- Map(unlist, mapply(extractList, total.ERs, dropList))

#'======================================================================================================================================
#' @description  Alternative solution for filtering overlap hit index given condition that we proposed

func <- function(peakset, ovHit, 
                 replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  require(purrr)
  stopifnot(length(peakset)>0)
  stopifnot(inherits(peakset[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  min.c <- ifelse(replicate.type=="Biological",
                  param <- length(peakset)-1,
                  param <- length(peakset))
  min.c <- as.integer(min.c)
  map(ovHit, lengths) %>% 
    reduce(`+`) %>% 
    map_lgl(`>=`, min.c) -> keep_me
  return(keep_me)
}

keep_me <- func(total.ERs, .hit_1, "Biological")
keepList <- map(.hit_1, ~.[keep_me])
dropList <- map(.hit_1, ~.[!keep_me])

##' @title .peakOverlapping
##' @description finding overlapped regions of current sample which supportedd with rest of replicate by pair-wise
##' @return Integer hit index
##' @example 

.peakOverlapping <- function(peakset, idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(idx))
  # set up the entry
  chosen <- peakset[[idx]]
  que.hit <- as(findOverlaps(chosen), "List")
  sup.hit <- lapply(peakset[- idx], function(ele_) {
    ans <- as(findOverlaps(chosen, ele_), "List")
    out.idx0 <- as(FUN(extractList(ele_$score, ans)), "List")
    out.idx0 <- out.idx0[!is.na(out.idx0)]
    ans <- ans[out.idx0]
  })
  res <- c(list(que.hit),sup.hit)
  names(res) <- c(names(peakset[idx]),names(peakset[-idx]))
  return(res)
}

##' @example
.hit_1 <- .peakOverlapping(peakset = total.ERs, idx = 1L, FUN = which.max)
.hit_2 <- .peakOverlapping(peakset = total.ERs, idx = 2L, FUN = which.max)
.hit_3 <- .peakOverlapping(peakset = total.ERs, idx = 3L, FUN = which.max)

##----------------------------------------------------------------------------------------------------------
#' @description make each hit table has same pattern

idx <- sort(names(.hit_1))
.hit_2 <- .hit_2[idx]
.hit_3 <- .hit_3[idx]

##----------------------------------------------------------------------------------------------------------
#' @description let hit table as matrix representation by using DataFrame
idx <- sort(names(.hit_1))
.hit_1 <- DataFrame(.hit_1)
.hit_2 <- DataFrame(.hit_2[idx])
.hit_3 <- DataFrame(.hit_3[idx])
Hit <- Map(rbind, .hit_1, .hit_2, .hit_3)

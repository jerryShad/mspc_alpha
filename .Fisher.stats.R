##' @title Fisher.stats
##' @description get combined pvalue for overlapped peaks by pair wise
##' @return GRanges
##' @example 

.Fisher.stats <- function(hitTB, peakset, verbose=FALSE, ...) {
  # input param checking
  stopifnot(inherits(peakset[[1]], "GRanges"))
  pval_List <- mapply(.get.pvalue, hitTB, peakset)
  .helper.PV <- function(p.list) {
    res <- sapply(p.list, function(x) {
      out <- ifelse(length(x)>0,
                    x,
                    0)
      out
    })
  }
  pval.TB <- Map(.helper.PV, pval_List)
  pval.TB <- data.frame(pval.TB)
  pv.stat <- rowSums( -2*log( pval.TB ), na.rm=T )
  Npresent <- rowSums( !is.na(pval.TB) )
  comb.pval <- pchisq( pv.stat, df=2*Npresent, lower=FALSE )
  res <- DataFrame(comb.pval)
  return(res)
}

.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

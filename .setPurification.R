#' @title .setPurification
#' @param peakList list of all confirmed and discarded peaks for each replicates
#' @param replicate.type indicate type of current sample
#' @return GRanges
#' @export
#' @importFrom dplyr setdiff
#' @author Julaiti Shayiding

.setPurification <- function(grs, replicate.type=c("Biological", "Technical"), ...) {
  # input param checking
  stopifnot(inherits(grs[[1]], "GRanges"))
  replicate.type = match.arg(replicate.type)
  DF <- lapply(grs, function(elm) {
    out <- as(elm, "data.frame")
  })
  if(replicate.type=="Biological") {
    res <- DF[[1]]
  } else {
    res <- setdiff(DF[[1]], DF[[2]])
  }
  out <- as(res, "GRanges")
  return(out)
}

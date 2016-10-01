# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling

## MSPC Package - Bioconductor Package for Multiple Sample Peak Calling
##
## Package workflow with all implemented function
##

##-----------------------------------------------------------------------------------
# 1 : read bed file as GRanges objects

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

# example
myData <- readPeakFiles(peakFolder = "data/")


##=========================================================================================
## 2: add pvalue as new metadata column for GRanges objects
.pvalueConversion <- function(x, pvalueBase = 1L, ...) {
  stopifnot(class(x) == "GRanges")
  stopifnot(is.numeric(pvalueBase))
  # explore score of all features
  if(is.null(x$pvalue)){
    x$p.value <- 10^(score(x)/(- pvalueBase))
    colnames(mcols(x))[3] <- "p.value"
  } else {
    x
  }
  return(x)
}

##==========================================================================================
## 3: remove all background noise from all replicates

.denoise_peakFiles <- function(peakFolder, denoise_threshold=1E-04, verbose=FALSE, ...) {
  if (verbose) {
    cat(">> filter out all background noise peaks whose pvalue above threshold \t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  stopifnot(is.numeric(denoise_threshold))
  if(!inherits(peakFolder[[1]], "GRanges")) {
    stop("file entry was not GRanges objects, invalid input")
  }
  res <- lapply(peakFolder, function(ele_) {
    ans <- .pvalueConversion(ele_, pvalueBase=1L)
    ans
  })
  ## TODO: FIXME & Beautify me as an elegant version
  filt <- lapply(res, function(ele_) {
    out <- subset(ele_, ele_$p.value < denoise_threshold)
    out
  })
  ## TO END;
  res <- filt
  return(res)
}

# testme:
all.peakFile <- .denoise_peakFiles(peakFolder = myData, denoise_threshold = 1e-4)

###===========================================================================================
### 4: peak overlapping all replicats simulatanously

.peakOverlapping <- function(peakset, idx=1L, FUN=which.min, ...) {
  # input param checking
  if(!inherits(peakset[[1]], "GRanges")) {
    stop("invalid input, type of entry must be GRanges objects")
  }
  stopifnot(is.numeric(idx))
  #FUN <- match.arg(FUN, which.min)
  
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

# testme:
# example (a.k.a, all needed test must be done):
test_1 <- .peakOverlapping(peakset = all.peakFile , idx = 1L, FUN = which.min)
test_2 <- .peakOverlapping(peakset = all.peakFile , idx = 2L, FUN = which.min)
test_3 <- .peakOverlapping(peakset = all.peakFile, idx = 3L, FUN = which.min)

#--------------------------------------------------------------------------------------
# TODO : implement one small utility function:

hit.list <- list(test1 = test_1, test2 = test_2, test3 = test_3)

####===================================================================================
#### 5: normalize all overlap hitlist into one


Hit2IntVec <- function(ovHit, verbose=FALSE, ...) {
  # input param cheching
  if(verbose) {
    cat(">> loading your hitlist...\t\t",
        format(Sys.time(), "%Y-%m-%d %X"), "\n")
  }
  res <- lapply(ovHit, function(ele_) {
    out <- drop(ele_)
    out[is.na(out)] <- 0L
    out
  })
  result <- do.call("cbind", res)
  return(result)
}

# testme:
hit.tmp <- Map(Hit2IntVec, hit.list)
.hitDF <- lapply(hitAsVect, as.data.frame)

library(dplyr)
hitComb <- bind_rows(.hitDF) %>% distinct
.hitMT <- as.matrix(hitComb)

.hitALL <- apply(.hitMT, 2, function(ele_) {
  ele_ <- as(ele_, "IntegerList")
  ele_[all(ele_==0L)] <- IntegerList(integer(0))
  ele_
})

#####=========================================================================================
#------ Brand new updated solution for filtering of sufficient overlap peak condition -------#

library(purrr)

min.c <- ifelse(replicate.type=="Biological",
                param <- length(peakset)-1,
                param <- length(peakset))

map(final_hit, lengths) %>%
  reduce(`+`) %>%
  map_lgl(`>=`, min.c) -> keep_me

keep_list <- map(final_hit, ~.[keep_me])
drop_list <- map(final_hit, ~.[!keep_me])

discarded_peaks <- mapply(extractList, peak.list, drop_list)
allDiscarded_0 <- lapply(discarded_peaks, function(ele_) {
  ans <- unlist(ele_)
  ans
})


###======================================================================================
.get.pvalue <- function(ovHit, obj, verbose=FALSE, ...) {
  # input param checking
  stopifnot(class(obj)=="GRanges")
  res <- extractList(obj$p.value, ovHit)
  return(res)
}

# testme:
pvl.list <- mapply(.get.pvalue, keep_list, peak.list)

#----------------------------------------------------------------------------------------
.helper.getPvalue <- function(pvlist, ...) {
  # input param checking
  res <- sapply(pvlist, function(x) {
    out <- ifelse(length(x)>0,
                  x,
                  0)
    out
  })
  return(res)
}

# testme:
list.pval <- data.frame(mapply(.helper.getPvalue, pvl.list))

#-----------------------------------------------------------------------------------------
# fisher method:
get.fisherScore <- function(pv.list, verbose=FALSE, ...) {
  # input param checking
  require(metap)
  Fisher.score <- suppressWarnings(
    out <- apply(pv.list[,], 1, function(ro) {
      ans <- sumlog(ro)$p
    })
  )
  res <- Fisher.score
  return(res)
}

comb.p <- get.fisherScore(list.pval)

#------------------------------------------------------------------------------------------
.helper.fishScore <- function(pvaList, verbose=FALSE, ...) {
  # input param checking
  fishScore <- get.fisherScore(pvaList)
  corr.vectDim <- lapply(pvaList, function(ele_) {
    tmp.df <- data.frame(cbind("p" = ele_, "comb.pvalue" = fishScore))
    out <- tmp.df[tmp.df$p != 0.000000e+00, ]
    out$p <- NULL
    out
  })
  res <- corr.vectDim
  return(res)
}
##
comb.pvl <- .helper.fishScore(pvaList = list.pval)

#-------------------------------------------------------------------------------------------

.hitAsGR <- function(gr, hit, verbose=FALSE, ...) {
  #input parama checking
  stopifnot(class(gr)=="GRanges")
  stopifnot(class(hit)=="CompressedIntegerList")
  .asGR <- gr[unlist(extractList(seq_along(gr), hit))]
  res <- .asGR
  return(res)
}

#--------------------------------------------------------------------------------------------
hit.peaks <- mapply(.hitAsGR, peak.list, keep_list)

candidate.peaks <- Map(cbind, hit.peaks, comb.pvl)

#---------------------------------------------------------------------------------------------

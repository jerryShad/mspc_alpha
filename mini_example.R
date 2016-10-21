# MSPC Project - Bioconductor Package for Multiple Sample Peak Calling

foo <- GRanges( seqnames=Rle("chr1", 4),ranges=IRanges(c(5,21,45,73), c(16,27,58,80)),
                rangeName=c("a1", "a2", "a3", "a4"), score=c(22, 6,13, 7))

bar <- GRanges(seqnames=Rle("chr1", 7),ranges=IRanges(c(7,14,26,39,49,65,77), c(10,19,34,43,53,71,91)),
               rangeName=c("b1", "b2", "b3", "b4", "b5", "b6","b7"), score=c(5, 9, 12, 5, 7, 7,16))

cat <- GRanges(seqnames=Rle("chr1", 5),ranges=IRanges(c(3,12,41,51,68), c(8,29,47,56,81)),
               rangeName=c("c1", "c2", "c3", "c4","c5"), score= c(4, 21, 7, 6, 15))

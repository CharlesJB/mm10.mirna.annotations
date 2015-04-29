#' Convert a GRanges from mm9 to mm10
#' 
#' @param gr The \code{GRanges} to convert.
#'
#' @return A \code{GRanges} object corresponding to the mm10 positions.
#'
#' @examples
#' \dontrun{
#'   gr_mm10 <- liftOver_mm9ToMm10(gr_mm9)
#' }
liftOver_mm9ToMm10 <- function(gr) {
  ch <- system.file("extdata/mm9ToMm10.over.chain",
                    package = "mm10.mirna.annotations")
  ch <- rtracklayer::import.chain(ch)
  GenomeInfoDb::seqlevelsStyle(gr) <- "UCSC"
  gr <- GenomeInfoDb::keepStandardChromosomes(gr)
  gr <- GenomicRanges::unlist(rtracklayer::liftOver(gr, ch))
  gr <- GenomeInfoDb::sortSeqlevels(gr)
  seqinfo <- GenomeInfoDb::Seqinfo(genome = "mm10")
  seqinfo <- suppressWarnings(GenomeInfoDb::keepStandardChromosomes(seqinfo))
  seqinfo <- GenomeInfoDb::sortSeqlevels(seqinfo)
  GenomeInfoDb::seqinfo(gr) <- seqinfo
  gr
}

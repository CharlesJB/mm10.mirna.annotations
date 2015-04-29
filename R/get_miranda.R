#' Prepare miRanda object
#' 
#' This function will download the mm10 miRanda database and convert it into
#' \code{GRanges} format.
#'
#' Requires an internet connection unless a filename corresponding to the
#' downloaded database is provided.
#'
#' By default, the "Good mirSVR score, Conserved miRNA" version of the database
#' is downloaded from the "August 2010 Release Downloads" section. See:
#' http://www.microrna.org/microrna/getDownloads.do
#'
#' @param filename The path to the miRanda download database in .txt format.
#'                 If \code{NULL}, the database will be downloaded.
#'                 Default: \code{NULL}
#'
#' @return A \code{GRanges} object corresponding to the mm10 miRanda database.
#'
#' @examples
#' \dontrun{
#'   miranda_gr <- miranda()
#' }
#'
#' @import GenomicRanges
#' @export
get_miranda <- function(filename = NULL) {
  # Check argument
  if (!is.null(filename)) {
    msg <- "Error in miranda("
    msg <- paste0(msg, as.character(filename))
    msg <- paste0(msg, ") : ")
    msg <- paste0(msg, "invalid 'filename' argument.\n")
    if (!is.character(filename)) {
      stop(paste0(msg, "'filename' must be a character string"))
    } else if (!file.exists(filename)) {
      stop(paste0(msg, "'filename' does not exists"))
    }
  }

  # Download database
  if (is.null(filename)) {
    filename <- "mouse_predictions_S_C_aug2010.txt.gz"
    url <- paste("http://cbio.mskcc.org/microrna_data", filename, sep = "/")
    download.file(url, filename)
  }

  # Load database
  miranda <- read.table(gzfile(filename), header = TRUE, comment.char =  "",
                        stringsAsFactors = FALSE, sep = "\t")

  # Extract data.frame
  coord <- gsub("\\[|\\]", "", miranda[["genome_coordinates"]])
  coord <- data.frame(do.call("rbind", lapply(strsplit(coord, "[:-]"),
                                function(x) x[c(2:4)])))
  colnames(coord) <- c("seqnames", "start", "end")
  i <- colnames(miranda) %in% c("mirna_name", "gene_symbol")
  miranda <- cbind(coord, miranda[,i])

  # Clean data.frame
  miranda[["start"]] <- suppressWarnings(
                          as.numeric(as.character(miranda[["start"]])))
  miranda[["end"]] <- suppressWarnings(
                       as.numeric(as.character(miranda[["end"]])))
  miranda[["seqnames"]] <- as.character(miranda[["seqnames"]])
  miranda[["seqnames"]] <- paste0("chr", miranda[["seqnames"]])

  ## We need to remove some strange genome_coordinates that caused NA, i.e.:
  ## [mm9:9:119034019-119034036,119036106-119036109:+]
  i <- !is.na(miranda[["start"]]) & !is.na(miranda[["end"]])
  miranda <- miranda[i,]

  tmp <- miranda[["start"]]
  i <- miranda[["start"]] > miranda[["end"]]
  miranda[["start"]][i] <- miranda[["end"]][i]
  miranda[["end"]][i] <- tmp[i]
  rm(tmp)

 # Create GRanges object
  miranda <- GenomicRanges::makeGRangesFromDataFrame(miranda,
          seqinfo = GenomeInfoDb::Seqinfo(genome = "mm9"),
          keep.extra.columns = TRUE)
  miranda <- GenomeInfoDb::keepStandardChromosomes(miranda)

  # liftOver
  liftOver_mm9ToMm10(miranda)
}

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
miranda <- function(filename = NULL) {
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
  miranda[["start"]] <- as.numeric(miranda[["start"]])
  miranda[["end"]] <- as.numeric(miranda[["end"]])
  miranda[["seqnames"]] <- as.character(miranda[["seqnames"]])
  miranda[["seqnames"]] <- paste0("chr", miranda[["seqnames"]])

  tmp <- miranda[["start"]]
  i <- miranda[["start"]] > miranda[["end"]]
  miranda[["start"]][i] <- miranda[["end"]][i]
  miranda[["end"]][i] <- tmp[i]
  rm(tmp)
  i <- !grepl("random", miranda[["seqnames"]])
  miranda <- miranda[i,]

  # Create GRanges object
  gr <- GenomicRanges::makeGRangesFromDataFrame(miranda,
          seqinfo = GenomeInfoDb::Seqinfo(genome = "mm10"), keep.extra.columns = TRUE)
  GenomeInfoDb::keepStandardChromosomes(gr)
}

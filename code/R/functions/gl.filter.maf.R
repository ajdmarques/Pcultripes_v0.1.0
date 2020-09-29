#' Filter loci according to minor allele frequency
#'
#' This function removes loci from genomic data that have a minor allele
#' frequency above or below a desired threshold.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param threshold \code{numeric} threshold minimum frequency permitted.
#'   Defaults to 0.05.
#'
#' @param sense \code{character} sense for apply threshold. Defaults to
#'  code{">"} such that alleles with a frequency below the argument to
#'  \code{threshold} are omitted.
#'
#' @param verbose \code{integer} should information be printed during
#'   execution?
#'
#' @return \code{\link[adegenet]{genlight}}
gl.filter.maf <- function(x, threshold = 0.05, sense = "<", v = 0) {
  # assert that arguments are valid
  assertthat::assert_that(inherits(x, "genlight"),
                          assertthat::is.number(threshold),
                          isTRUE(threshold <= 1.0),
                          isTRUE(threshold >= 0.0),
                          assertthat::is.string(sense),
                          assertthat::is.number(v))
  # calculate maf values
  mafs <- dartR::gl.report.maf(x)$maf
  # find loci to drop
  loci.drop <- do.call(sense, list(mafs, threshold))
  # print information if required
  if (v > 0.5) {
    cat("Starting gl.filter.maf: Filtering on minor allele frequency\n")
    cat(paste0("  Initial no. of loci = ", length(x@loc.names), "\n"))
    cat(paste0("  Removing loci with MAF ", sense , " ", threshold, "\n"))
    cat(paste0("  No. of loci deleted = ", sum(loci.drop), "\n"))
    cat("Completed gl.filter.maf\n")
  }
  # if none return original x object
  if (sum(loci.drop) < 0.5)
    return(x)
  # otherwise remove loci
  x <- x[, !loci.drop]
  x@other$loc.metrics <- x@other$loc.metrics[!loci.drop, ]
  # return reuslt
  x
}

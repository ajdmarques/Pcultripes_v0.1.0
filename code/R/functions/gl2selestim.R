#' Save SelEstim data
#'
#' This function exports input data for SelEstim.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param path \code{character} file path.
#'
#' @return invisible \code{TRUE} indicating success.
gl2selestim <- function(x, path) {
  # assert valdi arguments
  assertthat::assert_that(inherits(x, "genlight"), assertthat::is.string(path))
  # calculate allele frequencies
  af <- adegenet::tab(adegenet::genind2genpop(dartR::gl2gi(x)),
                      NA.method = "zero")
  # reformat matrix for selestim
  af2 <- t(sapply(seq(1, ncol(af), 2), function(i) c(t(af[, c(i, i + 1)]))))
  # save data to disk
  cat(adegenet::nPop(x), "\n", nrow(af2), "\n", sep = "", file = path)
  write.table(af2, file = path, append = TRUE, quote = TRUE, sep = " ",
                 eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                 col.names = FALSE)
  # return success
  invisible(TRUE)
}

#' Remove duplicate IDs from DartSEQ dataset
#'
#' This function removes samples with duplicate identifiers from a SNP
#' genotype files from DartSEQ.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @details Samples with the greatest degree of missing are excluded.
#'
#' @return \code{\link[adegenet]{genlight}} object.
remove_duplicate_ids <- function(x) {
  # assert argument validity
  assertthat::assert_that(inherits(x, "genlight"))
  # find duplicated ids
  dup_pos <- unname(which(duplicated(x@ind.names)))
  dup_ids <- unname(x@ind.names[duplicated(x@ind.names)])
  # temporarily assign deduplicated ids
  old_ids <- x@ind.names
  x@ind.names[dup_pos] <- paste0(x@ind.names[dup_pos], "_duplicate")
  id_data <- tibble::tibble(old_id = old_ids, new_id = x@ind.names)
  # recalc statistics
  cr <- dartR::gl.report.callrate(x, method = "ind", plot = FALSE, v = 2)
  cr <- rowSums(is.na(as.matrix(x))) / adegenet::nLoc(x)
  # iterate over duplicated ids and remove the most missing duplicated ids
  for (i in seq_along(dup_ids)) {
    ## find out which id to keep
    curr_ids <- id_data$new_id[which(id_data$old_id == dup_ids[i])]
    curr_cr <- cr[curr_ids]
    curr_keep_id <- curr_ids[which.max(curr_cr)]
    curr_remove_ids <- setdiff(curr_ids, curr_keep_id)
    ## remove other ids
    x <- dartR::gl.drop.ind(x, curr_remove_ids)
    id_data <- id_data[!id_data$new_id %in% curr_remove_ids, , drop = FALSE]
    ## set keep id to original id value
    x@ind.names[which(x@ind.names == curr_keep_id)] <- dup_ids[i]
  }
  # recalculate stats
  x <- dartR::gl.recalc.metrics(x, v = 2)
  # return result
  x
}

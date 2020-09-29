# read data
snp_admix_data <- readRDS(snp_ol_admix_path)

# assign samples to populations, and exclude samples with population certainty
# below threshold
## add populations to the data used for the population clustering analysis
snp_admix_pop_data <- lapply(seq_along(snp_raw_pop_data), function(i) {
  ## generate pop ids
  pop_id <- factor(as.character(apply(as.matrix(snp_pop_prob_data[[i]]),
                                      1, which.max)))
  pop_prob <- apply(as.matrix(snp_pop_prob_data[[i]]), 1, max)
  retain_sample <- pop_prob > 0
  ## create new genlight object
  x <- snp_raw_pop_data[[i]]
  x@pop <- pop_id
  ## subset out individuals with prop below threshold if needed
  if (sum(!retain_sample) > 0)
    x <- dartR::gl.keep.ind(x, x@ind.names[retain_sample])
  ## return result
  x
})

## add populations to the data used for the outlier loci analyses
snp_admix_data <- lapply(seq_along(snp_admix_data), function(i) {
  ## generate pop ids
  pop_id <- factor(as.character(apply(as.matrix(snp_pop_prob_data[[i]]),
                                      1, which.max)))
  pop_prob <- apply(as.matrix(snp_pop_prob_data[[i]]), 1, max)
  retain_sample <- pop_prob > 0
  ## create new genlight object
  x <- snp_admix_data[[i]]
  x@pop <- pop_id
  ## subset out individuals with prop below threshold if needed
  if (sum(!retain_sample) > 0)
    x <- dartR::gl.keep.ind(x, x@ind.names[retain_sample])
  x
})

## add populations to the spatial metadata
for (i in seq_along(snp_admix_pop_data)) {
  ## set pop ids
  pop_data <- tibble::tibble(
    order_id = snp_admix_pop_data[[i]]@ind.names,
    lineage = as.numeric(as.character(snp_admix_pop_data[[i]]@pop)))
  snp_metadata[[i]] <- dplyr::left_join(snp_metadata[[i]], pop_data,
                                        by = "order_id")
  attr(snp_metadata[[i]], "sf_column") <- "geometry"
}

# save results
snp_admix_pop_path <- "data/intermediate/snp_admix_pop.rds"
saveRDS(snp_admix_pop_data, snp_admix_pop_path, compress = "xz")
saveRDS(snp_admix_data, snp_ol_admix_path, compress = "xz")

# cleanup
rm(snp_admix_pop_path)


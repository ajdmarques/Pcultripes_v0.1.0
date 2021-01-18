source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("06"))

# load data
snp_ol_data <- readRDS(snp_ol_path)

# create neutral snp data set
snp_neutral_data <- lapply(seq_along(snp_ol_data), function(i) {
  x <- snp_ol_data[[i]]
  #! insert line to ensure outlier have the same number of loci
  pcadapt_outlier_loci[[1]] <- pcadapt_outlier_loci[[1]][pcadapt_outlier_loci[[1]]$loc_name %in% bayescan_outlier_loci[[1]]$loc_name,]
  adaptive_loci <- which(selestim_outlier_loci[[i]]$outlier &
                         bayescan_outlier_loci[[i]]$outlier &
                         pcadapt_outlier_loci[[i]]$outlier)
  loci_keep <- setdiff(seq_along(x@loc.names), adaptive_loci)
  x <- x[, loci_keep, drop = FALSE]
  x@other$loc.metrics <- x@other$loc.metrics[loci_keep, , drop = FALSE]
  x
})

# create adaptive snp data set
snp_adaptive_data <- lapply(seq_along(snp_ol_data), function(i) {
  x <- snp_ol_data[[i]]
  #! insert line to ensure outlier have the same number of loci
  pcadapt_outlier_loci[[1]] <- pcadapt_outlier_loci[[1]][pcadapt_outlier_loci[[1]]$loc_name %in% bayescan_outlier_loci[[1]]$loc_name,]
  adaptive_loci <- which(selestim_outlier_loci[[i]]$outlier &
                         bayescan_outlier_loci[[i]]$outlier &
                         pcadapt_outlier_loci[[i]]$outlier)
  cat("  number outlier loci = ", length(adaptive_loci), "\n")
  if (length(adaptive_loci) == 0) {
    x <- NULL
  } else {
    x <- x[, adaptive_loci, drop = FALSE]
    x@other$loc.metrics <- x@other$loc.metrics[adaptive_loci, , drop = FALSE]
  }
  x
})

# save data
snp_neutral_path <- "data/intermediate/snp_adaptive.rds"
saveRDS(snp_neutral_data, snp_neutral_path, compress = "xz")
snp_adaptive_path <- "data/intermediate/snp_neutral.rds"
saveRDS(snp_adaptive_data, snp_adaptive_path, compress = "xz")

# clean up
rm(snp_ol_data, snp_neutral_data, snp_adaptive_data)

# save session
session::save.session(session_path("07"), compress = "xz")

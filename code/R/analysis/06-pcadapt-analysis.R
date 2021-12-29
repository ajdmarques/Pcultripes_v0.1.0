source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("05"))

# load parameters
pcadapt_parameters <- "code/parameters/pcadapt.toml" %>%
                      RcppTOML::parseTOML() %>%
                     `[[`(MODE)

# load functions
source("code/R/functions/gl2bed.R")

# read data
snp_ol_admix_data <- readRDS(snp_ol_admix_path)

# create directory to save pcadapt results
unlink("data/intermediate/pcadapt-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/pcadapt-results", showWarnings = FALSE,
           recursive = TRUE)

# create input data for pcadapt
pcadapt_input_data <- lapply(seq_along(snp_ol_admix_data), function(i) {
  ## create directory to save pcadapt results
  curr_dir <- paste0("data/intermediate/pcadapt-results/",
                     spp_parameters$species_atlas_names[i])
  dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)
  ## save bed file with snp data
  gl2bed(snp_ol_admix_data[[i]], paste0(curr_dir, "/data.bed"))
  ## read bed file for pcadapt
  pcadapt::read.pcadapt(paste0(curr_dir, "/data.bed"), type = "bed")
})

# run initial analysis to determine suitable K
pcadapt_inital_run <- lapply(pcadapt_input_data, pcadapt::pcadapt,
                             K = pcadapt_parameters$initial_k,
                             method = pcadapt_parameters$method,
                             min.maf = pcadapt_parameters$min_maf,
                             pca.only = TRUE)

# determine suitable k for inference run
pcadapt_best_k <- lapply(pcadapt_inital_run, function(x) {
  nFactors::nScree(x$d)$Components$noc
})

# run pcadapt using suitable k
pcadapt_results <- lapply(seq_along(pcadapt_input_data), function(i) {
  pcadapt::pcadapt(pcadapt_input_data[[i]], K = pcadapt_best_k[[i]],
                   method = pcadapt_parameters$method,
                   min.maf = pcadapt_parameters$min_maf)
})

# write Mahalanobis distances to table for use in RDA
pcadapt_maha <- as.data.frame(pcadapt_results[[1]]$scores)
write.csv(pcadapt_maha, paste0(curr_dir,"/pcadapt_Mahalanobis.csv"))

# compile tables with loci names indicating if they are outlier or not
# according to the pcadapt analysis
pcadapt_outlier_loci <- lapply(seq_along(pcadapt_results), function(i) {
  tibble::tibble(loc_name = snp_ol_admix_data[[i]]@loc.names,
                 qvalue = qvalue::qvalue(pcadapt_results[[i]]$pvalues,
                                         pi0.method = "bootstrap")$qvalues,
                 outlier = qvalue < pcadapt_parameters$fdr) %>%
  dplyr::mutate(outlier = dplyr::if_else(is.na(outlier), FALSE, .$outlier))
})

# cleanup
rm(snp_ol_admix_data, pcadapt_best_k, pcadapt_inital_run, pcadapt_input_data, pcadapt_results)

# save session
session::save.session(session_path("06"), compress = "xz")

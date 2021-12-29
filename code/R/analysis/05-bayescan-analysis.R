source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("04"))

# load parameters
bayescan_parameters <- "code/parameters/bayescan.toml" %>%
                        RcppTOML::parseTOML() %>%
                       `[[`(MODE)

# load functions
source("code/R/functions/gl2bayescan.R")

# read data
snp_ol_data <- readRDS(snp_ol_path)

# create directory to save BayeScan results
unlink("data/intermediate/bayescan-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/bayescan-results", showWarnings = FALSE,
           recursive = TRUE)

# run BayeScan
bayescan_outputs <- lapply(seq_along(snp_ol_data), function(i) {
  ## create directory to save BayeScan results
  curr_dir <- paste0("data/intermediate/bayescan-results/",
                     spp_parameters$species_atlas_names[i])
  dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)
  ## save data in BayeScan format
  gl2bayescan(snp_ol_data[[i]], paste0(curr_dir, "/data.txt"))
  ## run BayeScan per number of replicates
  bayescan_path <- normalizePath(paste0("code/bayescan/BayeScan2.1/binaries",
                                        "/BayeScan2.1_linux64bits"))
  sapply(seq_len(bayescan_parameters$reps), function(j) {
    file.copy(paste0(curr_dir, "/data.txt"),
              paste0(curr_dir, "/data_run_", j, ".txt"))
    system(paste0(
      bayescan_path, " ",
      normalizePath(paste0(curr_dir, "/data_run_", j, ".txt")), " ",
      "-od ", curr_dir, " ",
      "-threads ", general_parameters$number_threads, " ",
      "-n ", sprintf("%i", bayescan_parameters$n), " ",
      "-thin ", sprintf("%i", bayescan_parameters$thin), " ",
      "-nbp ", sprintf("%i", bayescan_parameters$nbp), " ",
      "-pilot ", sprintf("%i", bayescan_parameters$pilot), " ",
      "-burn ", sprintf("%i", bayescan_parameters$burn),
      "-snp"))
    ## read BayeScan outputs
    paste0(curr_dir, "/data_run_", j, "_fst.txt")
  })
})

# extract mcmc diagnostic data
bayescan_trace_data <- lapply(bayescan_outputs, function(x) {
  curr_paths <- gsub("_fst.txt", ".sel", x, fixed = TRUE)
  suppressWarnings(
    plyr::ldply(seq_along(curr_paths), function(i)
      data.table::fread(curr_paths[i], data.table = FALSE) %>%
      dplyr::mutate(iteration = V1) %>%
      dplyr::mutate(run = i))) %>%
  dplyr::select(run, dplyr::everything()) %>%
  tibble::as_tibble()
})

# calculate the Gelman Rubin statistics using replicate runs
bayescan_gr_data <- lapply(bayescan_trace_data, function(x) {
  ## create mcmc coda object
  l <- lapply(split(x, x$run), function(y) {
    coda::mcmc(y$logL, start = min(y$iteration), end = max(y$iteration),
               thin = diff(y$iteration[seq_len(2)]))
  })
  l <- coda::mcmc.list(l)
  ## save plot with diagnostic
  png(paste0("data/intermediate/bayescan-results/",
             spp_parameters$species_atlas_names[i],
             "/traceplot.png"))
  coda::traceplot(l, smooth = FALSE)
  dev.off()
  ## calculate GR statistic
  gr <- coda::gelman.diag(l, autoburnin = FALSE)
  tibble::tibble(est = gr$psrf[[1]], upper_ci =  gr$psrf[[2]])
})

# import BayeScan results
bayescan_results <- lapply(seq_along(bayescan_outputs), function(i) {
  suppressWarnings(plyr::ldply(seq_along(bayescan_outputs[[i]]), function(j)
    data.table::fread(bayescan_outputs[[i]][[j]], data.table = FALSE,
                      col.names = c("V1", "prob", "log10(PO)", "qval",
                                    "alpha", "fst")) %>%
    dplyr::mutate(run = i,
                  loc_name = as.character(snp_ol_data[[i]]@loc.names)))) %>%
  dplyr::select(loc_name, dplyr::everything(), -V1) %>%
  setNames(c("loc_name", "prob", "log10_PO", "qval", "alpha", "fst", "run")) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(outlier = qval < bayescan_parameters$fdr)
})

# compile tables with loci names indicating if they are outlier or not
# according to the BayeScan analysis
bayescan_outlier_loci <- lapply(bayescan_results, function(x) {
  x %>%
  dplyr::group_by(loc_name) %>%
  dplyr::summarize(outlier = all(outlier), qvalue = max(qval)) %>%
  dplyr::ungroup()
})

# cleanup
rm(snp_ol_data)

# save session
session::save.session(session_path("05"), compress = "xz")

source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("03"))

# load parameters
selestim_parameters <- "code/parameters/selestim.toml" %>%
                        RcppTOML::parseTOML() %>%
                        `[[`(MODE)

# load functions
source("code/R/functions/gl2selestim.R")

# read data
snp_ol_data <- readRDS(snp_ol_path)

# create directory to save selestim plots
unlink("data/intermediate/selestim-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/selestim-results", showWarnings = FALSE,
           recursive = TRUE)

# run selestim
selestim_results <- lapply(seq_along(snp_ol_data), function(i) {
  ## create directory to save SelEstim results
  curr_dir <- paste0("data/intermediate/selestim-results/",
                     spp_parameters$species_atlas_names[i])
  dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)
  ## save data in SelEstim format
  gl2selestim(snp_ol_data[[i]], paste0(curr_dir, "/data.txt"))
  ## run SelEstim
  selestim_path <- normalizePath(paste0("code/selestim/",
                                        "SelEstim_1.1.7/src/selestim"))
  plyr::ldply(seq_len(selestim_parameters$reps), function(j) {
    # main run
    o <- system(paste0(
      selestim_path, " ",
      "-file ", normalizePath(paste0(curr_dir, "/data.txt")), " ",
      "-outputs ", paste0(curr_dir, "/run_", j), " ",
      "-threads ", general_parameters$number_threads, " ",
      "-thin ", sprintf("%i", selestim_parameters$thin), " ",
      "-npilot ", sprintf("%i", selestim_parameters$npilot), " ",
      "-lpilot ", sprintf("%i", selestim_parameters$lpilot), " ",
      "-burnin ", sprintf("%i", selestim_parameters$burnin), " ",
      "-seed ", sprintf("%i", sample.int(1e+5, 1)), " ",
      "-length ", sprintf("%i", selestim_parameters$length), " ",
      "-calibration >/dev/null"))
    # read in results and apply thresholds
    cal_data <- data.table::fread(paste0(curr_dir, "/run_", j,
                                         "/calibration/summary_delta.out"),
                                  data.table = FALSE)
    out_data <- data.table::fread(paste0(curr_dir, "/run_", j,
                                         "/summary_delta.out"),
                                  data.table = FALSE)
    threshold <- quantile(cal_data$KLD, (1 - selestim_parameters$limit))
    out_data$outlier <- out_data$KLD > threshold
    out_data$replicate <- j
    out_data$name <- snp_ol_data[[i]]@loc.names
    # return data
    out_data
  })
})

# read in diagnostic data
selestim_trace_data <- lapply(seq_along(snp_ol_data), function(i) {
  ## specify directory
  d <- paste0("data/intermediate/selestim-results/",
              spp_parameters$species_atlas_names[i])
  ## find files
  curr_paths <- paste0("data/intermediate/selestim-results/",
                        spp_parameters$species_atlas_names[i],
                       "/run_", seq_len(selestim_parameters$reps),
                       "/trace_lambda.out")
  suppressWarnings(
    plyr::ldply(seq_along(curr_paths), function(i)
      data.table::fread(curr_paths[i], data.table = FALSE) %>%
      dplyr::mutate(run = i))) %>%
  dplyr::select(run, dplyr::everything()) %>%
  tibble::as_tibble()
})

# calculate the Gelman Rubin statistics using replicate runs
selestim_gr_data <- lapply(selestim_trace_data, function(x) {
  ## create mcmc coda object
  l <- lapply(split(x, x$run), function(y) {
    coda::mcmc(y$lambda, start = min(y$iteration), end = max(y$iteration),
               thin = diff(y$iteration[seq_len(2)]))
  })
  l <- coda::mcmc.list(l)
  ## save plot with diagnostic
  png(paste0("data/intermediate/selestim-results/",
              spp_parameters$species_atlas_names[i],
             "/traceplot.png"))
  coda::traceplot(l, smooth = FALSE)
  dev.off()
  ## calculate GR statistic
  gr <- coda::gelman.diag(l, autoburnin = FALSE)
  tibble::tibble(est = gr$psrf[[1]], upper_ci =  gr$psrf[[2]])
})

# compile tables with loci names indicating if they are outlier or not
# according to the outflank analysis
selestim_outlier_loci <- lapply(selestim_results, function(x) {
  x %>%
  dplyr::rename(loc_name = name) %>%
  dplyr::group_by(loc_name) %>%
  dplyr::summarize(outlier = all(outlier), kld = mean(KLD)) %>%
  dplyr::ungroup()
})

# cleanup
rm(snp_ol_data)

# save session
session::save.session(session_path("04"), compress = "xz")

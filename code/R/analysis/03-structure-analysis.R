source('code/R/functions/session_path.R')
### Run the following in bash along with tmux seession
# R CMD BATCH --no-restore --no-save code/R/analysis/03-structure-analysis.R
#
### tmux notes
# tmux
# tmux list-sessions
# tmux attach -t (No.)
# Ctrl+B -> D # Detach
 
# restore session
try(session::restore.session(session_path("02")))

# load parameters
structure_parameters <- "code/parameters/structure.toml" %>%
                              RcppTOML::parseTOML() %>%
                             `[[`(MODE)
clumpp_parameters <- "code/parameters/clumpp.toml" %>%
                     RcppTOML::parseTOML() %>%
                    `[[`(MODE)

# load functions
source("code/R/functions/structure_analysis.R")

# read data
snp_raw_pop_data <- readRDS(snp_raw_pop_path)
snp_ol_data <- readRDS(snp_ol_path)

# run structure
#unlink("data/intermediate/structure-results", recursive = TRUE, force = TRUE)
#structure_outputs <- lapply(seq_along(snp_raw_pop_data), function(i)
  do.call(structure_analysis,
          list(
            x = snp_raw_pop_data[[i]],
            dir = paste0("data/intermediate/structure-results/",
                          spp_parameters$species_atlas_names[i]),
            numruns = structure_parameters$numruns,
            burnin = structure_parameters$burnin,
            numreps = structure_parameters$numreps,
            noadmix = structure_parameters$noadmix,
            admburnin = structure_parameters$admburnin,
            updatefreq = structure_parameters$updatefreq,
            freqscorr = structure_parameters$freqscorr,
            maxpops = structure_parameters$maxpops,
            threads = general_parameters$number_threads)))

# extract mcmc diagnostic data
structure_trace_data <- lapply(structure_outputs, function(x) {
  # replace output files with log files
  x <- gsub(".txt_f", ".log", x, fixed = TRUE)
  # split file paths by population number
  log_list <- split(x, dirname(x))
  # extract mcmc data
  plyr::ldply(seq_along(log_list), function(i) {
    plyr::ldply(seq_along(log_list[[i]]), function(j) {
      curr_file <- readLines(log_list[[i]][[j]])
      mcmc_matrix <- curr_file[seq(grep("Rep#", curr_file, fixed = TRUE)[1],
                                   grep("MCMC completed", curr_file,
                                        fixed = TRUE)[1])]
      mcmc_matrix <- mcmc_matrix[!grepl("Alpha", mcmc_matrix, fixed = TRUE)]
      mcmc_matrix <- mcmc_matrix[!grepl("BURNIN", mcmc_matrix, fixed = TRUE)]
      mcmc_matrix <- mcmc_matrix[!grepl("Burnin", mcmc_matrix, fixed = TRUE)]
      mcmc_matrix <- mcmc_matrix[!grepl("completed", mcmc_matrix, fixed = TRUE)]
      mcmc_matrix <- mcmc_matrix[!grepl("Rep", mcmc_matrix, fixed = TRUE)]
      mcmc_matrix <- mcmc_matrix[which(nchar(mcmc_matrix) > 0)]
      mcmc_matrix <- gsub(":", " ", mcmc_matrix, fixed = TRUE)
      mcmc_matrix <- data.table::fread(paste(mcmc_matrix, collapse = "\n"),
                                       sep = " ", header = FALSE,
                                       data.table = FALSE)
      # parse mcmc matrix column names
      mcmc_colnames <- grep("Est Ln", curr_file, fixed = TRUE, value = TRUE)[1]
      mcmc_colnames <- gsub("#", "", mcmc_colnames, fixed = TRUE)
      mcmc_colnames <- gsub(":", "", mcmc_colnames, fixed = TRUE)
      mcmc_colnames <- gsub("(", "", mcmc_colnames, fixed = TRUE)
      mcmc_colnames <- gsub(")", "", mcmc_colnames, fixed = TRUE)
      mcmc_colnames <- strsplit(mcmc_colnames, "  ")
      mcmc_colnames <- sapply(mcmc_colnames, trimws)
      mcmc_colnames <- mcmc_colnames[which(nchar(mcmc_colnames) > 0)]
      mcmc_colnames <- gsub(" ", ".", mcmc_colnames, fixed = TRUE)
      names(mcmc_matrix) <- mcmc_colnames
      # add in run number column
      curr_run <- gsub("run_", "", basename(log_list[[i]][[j]]), fixed = TRUE)
      curr_run <- gsub("k_", "", curr_run, fixed = TRUE)
      curr_run <- gsub(".log", "", curr_run, fixed = TRUE)
      curr_run <- strsplit(curr_run, "_", fixed = TRUE)[[1]][[2]]
      mcmc_matrix$replicate <- as.numeric(curr_run)
      # add in number of populations
      mcmc_matrix$k <- as.numeric(gsub("k_", "",
                                       basename(dirname(log_list[[i]][[j]])),
                                       fixed = TRUE))
      # return output
      mcmc_matrix
    })
  }) %>% tibble::as_tibble()
})

# calculate the Gelman Rubin statistics using replicate runs
structure_gr_data <- lapply(structure_trace_data, function(x) {
  plyr::ddply(x, "k", function(y) {
    l <- plyr::dlply(y, "replicate", function(z) {
      b <- structure_parameters$admburnin + structure_parameters$burnin
      coda::mcmc(z$Ln.Like, start = b + 1,
                 thin = structure_parameters$updatefreq)
    })
    gr <- coda::gelman.diag(coda::mcmc.list(l), autoburnin = FALSE)
    data.frame(est = gr$psrf[[1]], upper_ci =  gr$psrf[[2]])
  })  %>% tibble::as_tibble()
})

# calculate Evanno's K statistics per population number
structure_evannok_data <- lapply(structure_outputs, function(x) {
  # set filename to save Evanno K plot
  p <- paste0(normalizePath(dirname(dirname(x[[1]]))), "/evannok")
  # calculate statistics
  x %>%
  pophelper::readQ() %>%
  pophelper::tabulateQ() %>%
  pophelper::summariseQ() %>%
  {withr::with_package("ggplot2",
    pophelper::evannoMethodStructure(., exportplot = TRUE,
                                     outputfilename = p))} %>%
  tibble::as_tibble()
})

# run clump for the best supporting number of populations
#unlink("data/intermediate/clumpp-results", recursive = TRUE, force = TRUE)
#snp_pop_prob_data <- lapply(seq_along(structure_outputs), function(i) {
  # find best supported number of populations
  best_k <- structure_evannok_data[[i]]$k[which.max(
              structure_evannok_data[[i]]$deltaK)]
  # set folder to run clumpp analysis
  curr_dir <- paste0("data/intermediate/clumpp-results/",
                      spp_parameters$species_atlas_names[i])
  dir.create(curr_dir, showWarnings = FALSE, recursive = TRUE)
  # subset to only include files for best k
  curr_files <- structure_outputs[[i]][paste0("k_", best_k) ==
                                       basename(dirname(
                                         structure_outputs[[i]]))]
  # run clumpp analysis
  ## note that the code in pophelper does not appear to actually run
  ## clumpp, it merely sets up the directory, with param file and input data
  result <-
    curr_files %>%
    pophelper::readQ() %>%
    {withr::with_dir(curr_dir,
      pophelper::clumppExport(., parammode = clumpp_parameters$mode,
                              paramrep = clumpp_parameters$repeats))}
  ## find out where pophelper stored the data
  curr_clumpp_dir <- dirname(dir(curr_dir, "paramfile", recursive = TRUE,
                                 full.names = TRUE))
  ## manually read in the clumpp config file and set the W and S parameters
  curr_clump_config <- readLines(paste0(curr_clumpp_dir, "/paramfile"))
  curr_clump_config[startsWith("W ", curr_clump_config)] <-
    paste0("W ", clumpp_parameters$W)
  curr_clump_config[startsWith("S ", curr_clump_config)] <-
    paste0("S ", clumpp_parameters$S)
  curr_clump_config <- writeLines(curr_clump_config,
                                  paste0(curr_clumpp_dir, "/paramfile"))
  ## run clumpp
  file.copy("code/clumpp/CLUMPP_Linux64.1.1.2/CLUMPP",
            paste0(curr_clumpp_dir, "/CLUMPP_linux"))
  system(paste("chmod 777", paste0(curr_clumpp_dir, "/CLUMPP_linux")))
  withr::with_dir(curr_clumpp_dir, system("./CLUMPP_linux"))
  # read clumpp results
  out <- tibble::as_tibble(data.table::fread(dir(curr_clumpp_dir,
                                                 "^.*merged\\.txt",
                                                 full.names = TRUE),
                                             data.table = FALSE))
  out <- out[, seq(2, ncol(out) - 1), drop = FALSE]
  names(out) <- paste0("pop_", seq_len(ncol(out)))
  out
})

# assign samples to populations, and exclude samples with population certainty
# below threshold
## add populations to the data used for the population clustering analysis
snp_pop_data <- lapply(seq_along(snp_raw_pop_data), function(i) {
  ## generate pop ids
  pop_id <- factor(as.character(apply(as.matrix(snp_pop_prob_data[[i]]),
                                      1, which.max)))
  pop_prob <- apply(as.matrix(snp_pop_prob_data[[i]]), 1, max)
  retain_sample <- pop_prob > structure_parameters$probability_threshold
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
snp_ol_data <- lapply(seq_along(snp_ol_data), function(i) {
  ## generate pop ids
  pop_id <- factor(as.character(apply(as.matrix(snp_pop_prob_data[[i]]),
                                      1, which.max)))
  pop_prob <- apply(as.matrix(snp_pop_prob_data[[i]]), 1, max)
  retain_sample <- pop_prob > structure_parameters$probability_threshold
  ## create new genlight object
  x <- snp_ol_data[[i]]
  x@pop <- pop_id
  ## subset out individuals with prop below threshold if needed
  if (sum(!retain_sample) > 0)
    x <- dartR::gl.keep.ind(x, x@ind.names[retain_sample])
  x
})

## add populations to the spatial metadata
for (i in seq_along(snp_pop_data)) {
  ## set pop ids
  pop_data <- tibble::tibble(
    order_id = snp_pop_data[[i]]@ind.names,
    lineage = as.numeric(as.character(snp_pop_data[[i]]@pop)))
  snp_metadata[[i]] <- dplyr::left_join(snp_metadata[[i]], pop_data,
                                        by = "order_id")
  attr(snp_metadata[[i]], "sf_column") <- "geometry"
}

# save results
snp_pop_path <- "data/intermediate/snp_pop.rds"
saveRDS(snp_pop_data, snp_pop_path, compress = "xz")
saveRDS(snp_ol_data, snp_ol_path, compress = "xz")

source('code/R/analysis/03b-structure-analysis-admix.R')

# cleanup
rm(snp_pop_data, pop_data, snp_raw_pop_data)

# save session
session::save.session(session_path("03"), compress = "xz")


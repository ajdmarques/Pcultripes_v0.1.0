#' Structure Analysis
#'
#' Perform a Structure analysis.
#'
#' @param x \code{\link[adegenet]{genlight}} object.
#'
#' @param numruns \code{numeric} Number of replicate Structure runs. Defaults
#'   to 2.
#'
#' @param maxpops \code{numeric} Number of populations assumed. Defaults to
#'   a sequence ranging from 1 to 6.
#'
#' @param burnin \code{numeric} Length of burnin period. Defaults to 10000.
#'
#' @param numreps \code{numeric} Number of MCMC iterations for inference.
#'   Defaults to 20000.
#'
#' @param noadmix \code{logical} Do not use admixture model. Defaults to
#'   \code{FALSE}.
#'
#' @param admburnin \code{numeric} Initial period of burnin with admixture
#'   model. Defaults to 500.
#'
#' @param freqscorr \code{logical} Allele frequencies are correlated among
#'   populations? Defaults to \code{TRUE}.
#'
#' @param updatefreq \code{numeric} Frequency to store updates to log-likelihood
#'  for traceplots. Defaults to 200.
#'
#' @param dir \code{character} working directory. Defaults to a temporary
#'   directory.
#'
#' @param threads \code{integer} number of threads for computational processing.
#'   Defaults to 1.
#'
#' @param structure_path \code{character} path to structure executable file.
#'
#' @return
structure_analysis <- function(x, numruns = 2, maxpops = seq_len(6),
                               burnin = 10000, numreps = 20000,
                               noadmix = FALSE, admburnin = 500,
                               freqscorr = TRUE, updatefreq = 200,
                               dir = tempdir(), threads = 1,
                               structure_path = "code/structure/structure") {
  # Assert valid arguments
  assertthat::assert_that(inherits(x, "genlight"),
                          assertthat::is.count(numruns),
                          is.numeric(maxpops),
                          all(vapply(maxpops, assertthat::is.count,
                                     logical(1))),
                          assertthat::is.count(numreps),
                          assertthat::is.flag(noadmix),
                          assertthat::is.count(admburnin),
                          assertthat::is.flag(freqscorr),
                          assertthat::is.count(updatefreq),
                          assertthat::is.string(dir),
                          assertthat::is.count(threads),
                          file.exists(structure_path))
  # Initialization
  ## create directory
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(dir, "/input"), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0(dir, "/output"), showWarnings = FALSE, recursive = TRUE)
  sapply(paste0(dir, "/output/k_", maxpops), dir.create, showWarnings = FALSE,
         recursive = TRUE)
  # Save data to file
  ## save snp data
  dartR::gl2structure(x, outfile = paste0(dir, "/input/data.str"))
  # Save parameters to file
  ## main parameters
  mp <- readLines("code/templates/mainparams.txt")
  mp <- gsub("$$BURNIN$$", sprintf("%i", burnin), mp, fixed = TRUE)
  mp <- gsub("$$NUMREPS$$", sprintf("%i", numreps), mp, fixed = TRUE)
  writeLines(mp, paste0(dir, "/input/mainparams.txt"))
  ## extra parameters
  ep <- readLines("code/templates/extraparams.txt")
  ep <- gsub("$$NOADMIX$$", sprintf("%i", noadmix), ep, fixed = TRUE)
  ep <- gsub("$$FREQSCORR$$", sprintf("%i", freqscorr), ep, fixed = TRUE)
  ep <- gsub("$$UPDATEFREQ$$", sprintf("%i", updatefreq), ep, fixed = TRUE)
  ep <- gsub("$$ADMBURNIN$$", sprintf("%i", admburnin), ep, fixed = TRUE)
  writeLines(ep, paste0(dir, "/input/extraparams.txt"))
  # Run structure
  ## initialize cluster
  if (threads > 1) {
    cl <- parallel::makeCluster(threads, type = "FORK")
    doParallel::registerDoParallel(cl)
  }
  ## initialize data.frame to track structure runs/replicates
  md <- expand.grid(k = maxpops, replicate = seq_len(numruns))
  md$seed <- sample.int(1e+10, nrow(md))
  ## run analysis in parallel
  out <- plyr::llply(seq_len(nrow(md)), .parallel = threads > 1, function(i) {
    curr_output_dir <-
    cmd <- paste0(
      normalizePath(structure_path), " ",
       "-m \"", normalizePath(paste0(dir, "/input/mainparams.txt")), "\" ",
       "-e \"", normalizePath(paste0(dir, "/input/extraparams.txt")), "\" ",
       "-i \"", normalizePath(paste0(dir, "/input/data.str")), "\" ",
       "-K ", md$k[i], " ",
       "-L ", sprintf("%i", length(x@loc.names)), " ",
       "-N ", sprintf("%i", length(x@ind.names)), " ",
       "-D ", md$seed[i], " ",
       "-o \"", normalizePath(paste0(dir, "/output/k_", md$k[i])), "/k_",
          md$k[i], "_run_", md$replicate[i], ".txt\" ",
       "> \"", normalizePath(paste0(dir, "/output/k_", md$k[i])), "/k_",
          md$k[i], "_run_", md$replicate[i], ".log\" 2>&1")
    o <- system(cmd, intern = TRUE)
  })
  ## kill cluster
  if (threads > 1) {
    doParallel::stopImplicitCluster()
    cl <- parallel::stopCluster(cl)
  }
  # Exports
  ## return output file names
  dir(paste0(dir, "/output"), "^.*\\.txt_f$", recursive = TRUE, full.names = TRUE)
}

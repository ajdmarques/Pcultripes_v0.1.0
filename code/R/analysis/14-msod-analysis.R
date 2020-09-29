source('code/R/functions/session_path.R')
# restore session
session::restore.session(session_path("13"))

# load packages (because sf doesn't play nicely with dplyr otherwise...)
library(sf)
library(dplyr)
library(magrittr)

# load parameters
msod_parameters <- "code/parameters/msod.toml" %>%
  RcppTOML::parseTOML() %>%
  `[[`(MODE)

# load data
snp_admix_data <- readRDS(snp_ol_admix_path)
## Remove Valdemenco
x <- snp_metadata[[1]][grep('Valdemanco',snp_metadata[[1]]$locality_id),]$order_id
snp_admix_data[[1]] <- snp_admix_data[[1]][!snp_admix_data[[1]]$ind.names %in% x, ]
#snp_admix_data[[1]] <- snp_admix_data[[1]][snp_admix_data[[1]]$ind.names %in% snp_metadata[[1]]$order_id,]
# create directory to save rda results
unlink("data/intermediate/msod-results", force = TRUE, recursive = TRUE)
dir.create("data/intermediate/msod-results", showWarnings = FALSE,
           recursive = TRUE)

# run msod analysis
msod_results <- lapply(seq_along(snp_admix_data), function(i) {
  ## extract data
  curr_g_data <- snp_admix_data[[i]]
  curr_m_data <- snp_metadata[[i]]
  ## add small jitter (<= 1 m) to coordinates to avoid duplicates
  curr_coords <- as(curr_m_data, "Spatial")@coords
  curr_coords[, 1] <- curr_coords[, 1] + runif(nrow(curr_coords))
  curr_coords[, 2] <- curr_coords[, 2] + runif(nrow(curr_coords))
  ## calculate MEMs
  curr_nb <- spdep::graph2nb(spdep::gabrielneigh(curr_coords), sym = TRUE)
  curr_listw <- spdep::nb2listw(curr_nb, style = "W")
  curr_d <- spdep::nbdists(curr_nb, curr_coords)
  curr_d <- lapply(curr_d, function(x) 1 / x)
  curr_listw <- spdep::nb2listw(curr_nb, glist = curr_d, style = "W")
  curr_mem <- adespatial::scores.listw(curr_listw, MEM.autocor = "all")
  curr_mem <- list(vectors = as.matrix(curr_mem),
                   values = attr(curr_mem, "values"))
  curr_mem$values <- curr_mem$values / abs(sum(curr_mem$values))
  ## correlate MEMs with loci
  curr_cor <- cor(as.matrix(snp_admix_data[[i]]), curr_mem$vectors,
                  use = "pairwise.complete.obs")
  ## calculate scores
  curr_cor_sq <- curr_cor ^ 2
  curr_cor_sq_centered <- apply(curr_cor_sq, 2, mean, na.rm = T)
  curr_dev <- sweep(curr_cor_sq, 2, curr_cor_sq_centered, "/") - 1
  curr_dev[curr_dev > 0] <- 0
  curr_dev <- apply(curr_dev, 1, sum)
  curr_z <- (curr_dev - mean(curr_dev)) / sd(curr_dev)
  ## return absolute values of z-scores
  abs(curr_z)
})

# compile tables with loci names indicating if they are outlier or not
# according to the MSOD analysis
msod_outlier_loci <- lapply(seq_along(snp_admix_data), function(i) {
  tibble::tibble(loc_name = snp_admix_data[[i]]@loc.names,
                 z_score = msod_results[[i]],
                 outlier = z_score > msod_parameters$threshold)
})

# save results
msod_results_path <- "data/intermediate/msod_results.rds"
saveRDS(msod_results, msod_results_path, compress = "xz")

# clean up
rm(snp_admix_data, msod_results, bayescan_results, lfmm_results, load.rda, 
   p.values, pcadapt_input_data, pcadapt_inital_run, pcadapt_results, pvalues, sambada_results,
   sambada_snp, scaleY, selestim_results, snp_admix_pop_data, snp_genind, snp_inp, snp_lfmm,
   snp_ol_data, snp_raw_pop_data, test, y)

# unload packages
detach("package:dplyr", unload = TRUE)
detach("package:sf", unload = TRUE)

# save session
session::save.session(session_path("14"), compress = "xz")

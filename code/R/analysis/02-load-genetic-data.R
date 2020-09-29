source('code/R/functions/session_path.R')
source('code/R/functions/remove_duplicate_ids.R')
# restore session
session::restore.session(session_path("01"))

# load parameters
genetic_data_parameters <-
  RcppTOML::parseTOML("code/parameters/genetic-data.toml") %>%
  `[[`(MODE)

# load packages (because sf doesn't play nicely with dplyr otherwise...)
library(sf)
library(dplyr)

# load functions
source("code/R/functions/gl.filter.maf.R")
source("code/R/functions/read_sample_data.R")

# load sample location data
snp_metadata <-
  lapply(spp_parameters$species_snp_metadata_path, read_pc_sample_data) %>%
  lapply(function(x) {
  sf::st_as_sf(x, coords = c("longitude", "latitude"), crs = 4326,
               agr = "constant") %>%
  sf::st_transform(general_parameters$crs)
})

# load genetic data
snp_raw_data <- lapply(spp_parameters$species_snp_raw_path,
                       dartR::gl.read.dart)

snp_clean_data <-
  lapply(seq_along(snp_raw_data), function(i) {
    cat("preliminary data cleaning:",
        basename(spp_parameters$species_snp_raw_path[i]), "\n")
    # initialization
    x <- snp_raw_data[[i]]
    # manually replace spaces in the sample + loci names with underscores
    x@loc.names <- gsub(" ", "_", x@loc.names, fixed = TRUE)
    x@ind.names <- gsub(" ", "_", x@ind.names, fixed = TRUE)
    # add GVA prefix to ids that are not in the metadata
    p1 <- which(!x@ind.names %in% snp_metadata[[i]]$order_id)
    nn <- paste0("GVA", x@ind.names[p1])
    p2 <- nn %in% snp_metadata[[i]]$order_id
    if (any(nn %in% snp_metadata[[i]]$order_id))
      x@ind.names[p1[p2]] <- nn[p2]
    # add IMS_ prefix to ids that are not in the metadata
    p1 <- which(!x@ind.names %in% snp_metadata[[i]]$order_id)
    nn <- paste0("IMS_", x@ind.names[p1])
    p2 <- nn %in% snp_metadata[[i]]$order_id
    if (any(nn %in% snp_metadata[[i]]$order_id))
      x@ind.names[p1[p2]] <- nn[p2]
    # remove GVA prefix from ids that are not in the metadata
    p1 <- which(!x@ind.names %in% snp_metadata[[i]]$order_id)
    nn <- gsub("GVA", "", x@ind.names[p1], fixed = TRUE)
    p2 <- nn %in% snp_metadata[[i]]$order_id
    if (any(nn %in% snp_metadata[[i]]$order_id))
      x@ind.names[p1[p2]] <- nn[p2]
    # remove duplicated identifiers by removing samples with most missingness
#    if (i == 3) {
#      x@ind.names[which(x@ind.names == "45908")] <- "8932"
#            x <- remove_duplicate_ids(x)
#    } else if (i == 2) {
      x@ind.names[which(x@ind.names == "108541")] <- "106854"
      x <- remove_duplicate_ids(x)
#    }
    # manually replace sample names where they aren't in the metadata
#    p1 <- which(!x@ind.names %in% snp_metadata[[i]]$order_id)
#    if ((length(p1) > 0) &&
#        assertthat::has_name(snp_metadata[[i]], "collector_id")) {
#      x@ind.names[p1] <- snp_metadata[[i]]$order_id[
#        match(gsub("_", " ", x@ind.names[p1], fixed = TRUE),
#              snp_metadata[[i]]$collector_id)]
#   }
    # assert that all individuals in the SNP dataset are in the metadata table
#    assertthat::assert_that(all(x@ind.names %in% snp_metadata[[i]]$order_id),
#                            anyDuplicated(x@ind.names) == 0)
    # clean data
    x %>%
      dartR::gl.filter.repavg(t = genetic_data_parameters$repavg_threshold,
                              v = 3) %>%
      dartR::gl.filter.callrate(
        method = "loc",
        threshold = genetic_data_parameters$locus_callrate_threshold,
        v = 3) %>%
      dartR::gl.filter.secondaries(method =
                                     genetic_data_parameters$sec_method,
                                   v = 3) %>%
      dartR::gl.filter.maf(threshold = genetic_data_parameters$min_maf, v = 1) %>%
      dartR::gl.filter.callrate(
        method = "ind",
        threshold = genetic_data_parameters$sample_callrate_threshold,
        v = 3)
  })

# filter loci if in debugging mode to keep only first 1000 loci
if (identical(MODE, "debug")) {
  snp_clean_data <-
    lapply(seq_along(snp_clean_data), function(i) {
      ## extract data
      x <- snp_clean_data[[i]]
      all_n <- x@loc.names
      n_keep <- 5000
      ## keep only first bunch of loci if there's more than that
      if (length(all_n) > n_keep) {
        keep_n <- sample(all_n, min(n_keep, length(all_n)))
        drop_n <- setdiff(all_n, keep_n)
        x <- dartR::gl.drop.loc(x, drop_n)
      }
      ## return object
      x
    })
}

# filter data for pop analysis
snp_raw_pop_data <-
  lapply(seq_along(snp_clean_data), function(i) {
    cat("data prep for population analysis:",
        spp_parameters$species_sci_names[i], "\n")
    {snp_clean_data[[i]]} %>%
    `slot<-`("pop", value = factor(rep("all", length(.@ind.names)))) %>%
    dartR::gl.filter.hwe(alpha = genetic_data_parameters$hwe_alpha) %>%
    `slot<-`("pop", value = NULL)
})

# filter data for outlier loci analysis
snp_ol_data <-
  lapply(seq_along(snp_raw_data), function(i) {
    cat("data prep for outlier loci analysis:",
        spp_parameters$species_sci_names[i], "\n")
    {snp_clean_data[[i]]}
})

# assert that both population and outlier loci data sets have same samples
lapply(seq_along(snp_raw_data), function(i) {
  assertthat::assert_that(identical(snp_ol_data[[i]]@ind.names,
                                    snp_raw_pop_data[[i]]@ind.names))
})

# assert that each SNP dataset contains at least one SNP
lapply(seq_along(snp_raw_data), function(i) {
  assertthat::assert_that(length(snp_ol_data[[i]]@ind.names) > 0,
                          length(snp_raw_pop_data[[i]]@ind.names) > 0,
                          length(snp_ol_data[[i]]@loc.names) > 0,
                          length(snp_raw_pop_data[[i]]@loc.names) > 0)
})

# save data
snp_raw_path <- "data/intermediate/snp_raw.rds"
saveRDS(snp_raw_data, snp_raw_path, compress = "xz")
snp_clean_path <- "data/intermediate/snp_clean.rds"
saveRDS(snp_clean_data, snp_clean_path, compress = "xz")
snp_raw_pop_path <- "data/intermediate/snp_raw_pop.rds"
saveRDS(snp_raw_pop_data, snp_raw_pop_path, compress = "xz")
snp_ol_path <- "data/intermediate/snp_ol.rds"
saveRDS(snp_ol_data, snp_ol_path, compress = "xz")
snp_ol_admix_path <- "data/intermediate/snp_ol_admix.rds"   ## Aditional save for admixed data
saveRDS(snp_ol_data, snp_ol_admix_path, compress = "xz")            ## Aditional save for admixed data

# cleanup
rm(snp_raw_data, snp_ol_data, snp_raw_pop_data, snp_clean_data)
if (file.exists("Rplots.pdf"))
  unlink("Rplots.pdf")

# unload packages
detach("package:dplyr", unload = TRUE)
detach("package:sf", unload = TRUE)

# save session
session::save.session(session_path("02"), compress = "xz")


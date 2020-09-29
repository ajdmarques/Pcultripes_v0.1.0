#' Read Hyla molleri sample data
#'
#' This function reads and formats the sample data for \emph{Hyla molleri}.
#'
#' @param x \code{character} file path.
#'
#' @return \code{data.frame} object.
read_hm_sample_data <- function(x) {
  # assert arguments are valid
  assertthat::assert_that(assertthat::is.string(x), file.exists(x))
  # format data
  readxl::read_excel(x) %>%
    dplyr::select(1:12) %>%
    setNames(c("cluster_id", "locality_id", "order_id", "genus", "species_name",
               "sex", "province", "locality_description", "latitude",
               "longitude", "country", "collector_id")) %>%
    dplyr::filter(!is.na(order_id)) %>%
    dplyr::mutate(sex = dplyr::case_when(sex == "Hembra" ~ "female",
                                         sex == "Macho" ~ "male",
                                         TRUE ~ NA_character_)) %>%
    dplyr::mutate(cluster_id = zoo::na.locf(cluster_id)) %>%
    dplyr::mutate(locality_id = zoo::na.locf(locality_id)) %>%
    dplyr::mutate(locality_id = paste0(cluster_id, "_", locality_id)) %>%
    dplyr::mutate(sci_name = paste(genus, species_name)) %>%
    dplyr::select(-genus, -species_name)
}

#' Read Pelobates cultripes sample data
#'
#' This function reads and formats the sample data for
#' \emph{Pelobates cultripes}.
#'
#' @param x \code{character} file path.
#'
#' @return \code{data.frame} object.
read_pc_sample_data <- function(x) {
  # assert arguments are valid
  assertthat::assert_that(assertthat::is.string(x), file.exists(x))
  # format data
  readxl::read_excel(x, col_types = "text") %>%
    dplyr::select(1:4, 8) %>%
    setNames(c("order_id", "locality_id", "latitude", "longitude",
               "cluster_id")) %>%
    dplyr::mutate(latitude = as.numeric(latitude),
                  longitude = as.numeric(longitude)) %>%
    dplyr::select(order_id, locality_id, longitude, latitude, cluster_id)
}

#' Read Rana iberica sample data
#'
#' This function reads and formats the sample data for \emph{Rana iberica}.
#'
#' @param x \code{character} file path.
#'
#' @return \code{data.frame} object.
read_ri_sample_data <- function(x) {
  # assert arguments are valid
  assertthat::assert_that(assertthat::is.string(x), file.exists(x))
  # format data
  readxl::read_excel(x, col_types = "text", skip = 1) %>%
    dplyr::select(1, 2, 3, 4, 5) %>%
    setNames(c("order_id", "locality_id", "cluster_id",
               "latitude", "longitude")) %>%
    dplyr::mutate(
      order_id = gsub(" ", "_", order_id, fixed = TRUE),
      latitude = as.numeric(latitude),
      longitude = as.numeric(longitude)) %>%
    dplyr::select(order_id, locality_id, longitude, latitude, cluster_id)
}

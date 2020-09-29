source('code/R/functions/session_path.R')
# restore session
try(session::restore.session(session_path("00")))

# load packages (because sf doesn't play nicely with dplyr otherwise...)
library(sf)
library(dplyr)
library(magrittr)

# load parameters
spp_parameters <- RcppTOML::parseTOML("code/parameters/species.toml") %>%
                    `[[`(MODE)
assertthat::assert_that(
  length(spp_parameters$species_atlas_names) ==
    length(spp_parameters$species_sci_names),
  length(spp_parameters$species_atlas_names) ==
    length(spp_parameters$species_eng_names))

# create planning units
# load atlas data
atlas_data <- sf::read_sf("data/raw/atlas/Atlas.shp") %>%
              lwgeom::st_make_valid()

# prepare atlas data
## subset to species
atlas_data <- atlas_data[, spp_parameters$species_atlas_names, drop = FALSE]

## reproject to coordinate system
atlas_data <- atlas_data %>%
              sf::st_transform(general_parameters$crs)

# make study area raster from the polygon data for spatial conformity
## create empty raster
study_area_raster_data <-
  raster::raster(
    xmn = raster::xmin(raster::extent(atlas_data)),
    ymn = raster::ymin(raster::extent(atlas_data)),
    xmx = plyr::round_any(raster::xmax(raster::extent(atlas_data)),
                           general_parameters$resolution, ceiling),
    ymx = plyr::round_any(raster::ymax(raster::extent(atlas_data)),
                           general_parameters$resolution, ceiling),
    resolution = rep(general_parameters$resolution, 2),
    crs = sp::CRS(general_parameters$crs))
study_area_raster_data[] <- seq_len(raster::ncell(study_area_raster_data))

## find out which raster cells are inside grid cells
cells <- raster::as.data.frame(study_area_raster_data, na.rm = FALSE, xy = TRUE)
cells <- sp::SpatialPointsDataFrame(cells[, c("x", "y")], cells,
                           proj4string = sp::CRS(general_parameters$crs)) %>%
         as("sf")
cells <- sf::st_intersection(cells, sf::st_union(sf::st_geometry(atlas_data)))
study_area_raster_data[] <- dplyr::if_else(study_area_raster_data[] %in%
                                           cells$layer, 1, NA_real_)

# create RasterStack of species presence/absences
## initialize data
habitat_raster_data <-
  study_area_raster_data[[rep(1, length(spp_parameters$species_atlas_names))]]
names(habitat_raster_data) <- spp_parameters$species_sci_names
## add in presence absence data
for (i in seq_along(spp_parameters$species_sci_names)) {
  ## extract grid centroids for where species found
  pts <- {atlas_data[atlas_data[[i]] > 0.5, ]} %>%
         sf::st_centroid() %>%
         as("Spatial")
  ## extract cell indices
  cells <- raster::cellFromXY(study_area_raster_data, pts)
  cells <- cells[is.finite(study_area_raster_data[cells])]
  ## reset all nonNA cells to zero
  habitat_raster_data <- habitat_raster_data[[i]] * 0
  ## add in presences
  habitat_raster_data[cells] <- 1
}

# mask out planning units that are not occupied by any species
study_area_raster_data[sum(habitat_raster_data) == 0] <- NA_real_
habitat_raster_data %<>% raster::mask(study_area_raster_data)

# create raster with the proportion of land in each pixel
## load land data
## extract boundary for mainland of Iberian Peninsula
ip_data <- dir("data/raw/gadm", "^.*\\.rds$", full.names = TRUE) %>%
           lapply(readRDS) %>%
           do.call(what = rbind) %>%
           sf::st_transform(general_parameters$crs) %>%
           sf::st_union() %>%
           as("Spatial") %>%
           sp::disaggregate() %>%
           as("sf") %>%
           mutate(area = sf::st_area(.)) %>%
           arrange(dplyr::desc(area)) %>%
           filter(seq_along(area) == 1)

## create land raster
land_raster_data <-
  ip_data %>%
  fasterize::fasterize(raster::disaggregate(study_area_raster_data, 10),
                       background = 0) %>%
  raster::aggregate(fact = 10, fun = mean) %>%
  raster::mask(study_area_raster_data)

# specify cost data as percentage of land inside each planning unit
cost_raster_data <- land_raster_data

# mask out planning units with no data
study_area_raster_data %<>% raster::mask(sum(raster::stack(
  habitat_raster_data, cost_raster_data, land_raster_data), na.rm = FALSE))
habitat_raster_data %<>% raster::mask(study_area_raster_data)
cost_raster_data %<>% raster::mask(study_area_raster_data)
land_raster_data %<>% raster::mask(study_area_raster_data)

# create raster with planning unit identifiers (1--number of pixels)
planning_unit_raster_data <-
  study_area_raster_data
planning_unit_raster_data[!is.na(planning_unit_raster_data)] <-
 seq_len(raster::cellStats(planning_unit_raster_data, "sum"))

# verify that raster data conform to same spatial grid
assertthat::assert_that(
  raster::compareRaster(study_area_raster_data, planning_unit_raster_data,
                        res = TRUE, tolerance = 1e-5),
  identical(raster::Which(is.na(study_area_raster_data), cells = TRUE),
            raster::Which(is.na(planning_unit_raster_data), cells = TRUE)),
  identical(raster::Which(!is.na(study_area_raster_data), cells = TRUE),
            raster::Which(!is.na(planning_unit_raster_data), cells = TRUE)),
  raster::compareRaster(study_area_raster_data, habitat_raster_data[[1]],
                        res = TRUE, tolerance = 1e-5),
  identical(raster::Which(is.na(study_area_raster_data), cells = TRUE),
            raster::Which(is.na(habitat_raster_data[[1]]), cells = TRUE)),
  identical(raster::Which(!is.na(study_area_raster_data), cells = TRUE),
            raster::Which(!is.na(habitat_raster_data[[1]]), cells = TRUE)),
  raster::compareRaster(study_area_raster_data, cost_raster_data,
                        res = TRUE, tolerance = 1e-5),
  identical(raster::Which(is.na(study_area_raster_data), cells = TRUE),
            raster::Which(is.na(cost_raster_data), cells = TRUE)),
  identical(raster::Which(!is.na(study_area_raster_data), cells = TRUE),
            raster::Which(!is.na(cost_raster_data), cells = TRUE)))

# save data
ip_path <- "data/intermediate/ip.shp"
sf::st_write(ip_data, ip_path, delete_dsn = TRUE)
study_area_path <- "data/intermediate/study_area.tif"
raster::writeRaster(study_area_raster_data, study_area_path, NAflag = -9999,
                    overwrite = TRUE)
cost_path <- "data/intermediate/cost.tif"
raster::writeRaster(cost_raster_data, cost_path, NAflag = -9999,
                    overwrite = TRUE)
habitat_path <- "data/intermediate/habitat.tif"
raster::writeRaster(habitat_raster_data, habitat_path, NAflag = -9999,
                    overwrite = TRUE)
planning_unit_path <- "data/intermediate/pu.tif"
raster::writeRaster(planning_unit_raster_data, planning_unit_path,
                    NAflag = -9999, overwrite = TRUE)
land_path <- "data/intermediate/land.tif"
raster::writeRaster(land_raster_data, land_path,
                    NAflag = -9999, overwrite = TRUE)
# cleanup
rm(ip_data, study_area_raster_data, planning_unit_raster_data, cost_raster_data,
   habitat_raster_data, land_raster_data, atlas_data)

# unload packages
try(detach("package:dplyr", unload = TRUE))
try(detach("package:sf", unload = TRUE))

# save session
session::save.session(session_path("01"), compress = "xz")

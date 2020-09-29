# This script generates fake data for Rana and Pelobates by copying the Hyla
# data and randomizing which samples are associated with which genetic data

# Initialization
## load packages
library(magrittr)

## set paths
input_path <- paste0("/home/jeff/GitHub/amph-evol-conservation/data/raw/",
                     "hyla_dart2017/OrderAppendix_1_DHyl17-2984/",
                     "Report_DHyl17-2984_1_moreOrders_SNP_2.csv")

output_paths <- c(
  paste0("/home/jeff/GitHub/amph-evol-conservation/data/raw/",
         "pelobates_dart2019_fake/OrderAppendix_1_DHyl17-2984/",
         "Report_DHyl17-2984_1_moreOrders_SNP_2.csv"),
  paste0("/home/jeff/GitHub/amph-evol-conservation/data/raw/",
         "rana_dart2019_fake/OrderAppendix_1_DHyl17-2984/",
         "Report_DHyl17-2984_1_moreOrders_SNP_2.csv"))

# Preliminary processing
## load data
input_data <- readLines(input_path) %>%
              strsplit(",") %>%
              {matrix(unlist(.), byrow = TRUE, nrow = length(.),
                      ncol = length(.[[1]]))}

# Main processing
## find sample name cells
row <- min(which(input_data[, 1] != "*"))
cols <- seq(min(which(input_data[1, ] != "*")), ncol(input_data))

# Exports
lapply(output_paths, function(x) {
  z <- input_data
  z[row, cols] <- sample(z[row, cols])
  write.table(z, x, row.names = FALSE, col.names = FALSE, quote = FALSE,
              sep = ",")
  TRUE
})

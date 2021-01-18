# set seed for reproducibility
source('code/R/functions/session_path.R')
set.seed(500)

# define settings for analysis
MODE = "release"

# load packages
library(magrittr)

# load parameters
general_parameters <- RcppTOML::parseTOML("code/parameters/general.toml")
general_parameters <- general_parameters[[MODE]]

# load miscellaneous finstunctions
source("code/R/functions/misc.R")

# save session
session::save.session(session_path("00"), compress = "xz")
 
